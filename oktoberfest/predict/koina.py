import time
import warnings
from functools import partial
from typing import Dict, Generator, KeysView, List, Optional, Union

import numpy as np
import pandas as pd
from tqdm.auto import tqdm
from tritonclient.grpc import (
    InferenceServerClient,
    InferenceServerException,
    InferInput,
    InferRequestedOutput,
    InferResult,
)

alternative_column_map = {
    "peptide_sequences": "MODIFIED_SEQUENCE",
    "precursor_charges": "PRECURSOR_CHARGE",
    "collision_energies": "COLLISION_ENERGY",
    "fragmentation_types": "FRAGMENTATION",
    "instrument_types": "INSTRUMENT_TYPES",
}


class Koina:
    """A class for interacting with Koina models for inference."""

    model_inputs: Dict[str, str]
    model_outputs: Dict[str, np.ndarray]
    batch_size: int
    _response_dict: Dict[int, Union[InferResult, InferenceServerException]]

    def __init__(
        self,
        model_name: str,
        server_url: str = "koina.proteomicsdb.org:443",
        ssl: bool = True,
        targets: Optional[List[str]] = None,
        disable_progress_bar: bool = False,
    ):
        """
        Initialize a KoinaModel instance with the specified parameters.

        This constructor initializes the KoinaModel instance, connecting it to the specified Inference Server.
        It checks the availability of the server, the specified model, retrieves input and output information,
        and determines the maximum batch size supported by the model's configuration.
        Note: To use this class, ensure that the inference server is properly configured and running,
        and that the specified model is available on the server.

        :param model_name: The name of the Koina model to be used for inference.
        :param server_url: The URL of the inference server. Defaults to "koina.proteomicsdb.org:443".
        :param ssl: Indicates whether to use SSL for communication with the server. Defaults to True.
        :param targets: An optional list of targets to predict. If this is None, all model targets are
            predicted and received.
        :param disable_progress_bar: Whether to disable the progress bar showing the progress of predictions.
        """
        self.model_inputs = {}
        self.model_outputs = {}
        self._response_dict = {}

        self.model_name = model_name
        self.url = server_url
        self.ssl = ssl
        self.disable_progress_bar = disable_progress_bar
        self.client = InferenceServerClient(url=server_url, ssl=ssl)

        self.type_convert = {
            "FP32": np.dtype("float32"),
            "BYTES": np.dtype("O"),
            "INT16": np.dtype("int16"),
            "INT32": np.dtype("int32"),
            "INT64": np.dtype("int64"),
        }

        self._is_server_ready()
        self._is_model_ready()

        self.__get_inputs()
        self.__get_outputs(targets)
        self.__get_batchsize()

    @property
    def response_dict(self):
        """The dictionary containing raw InferenceResult/InferenceServerException objects (values) for a given request_id (key)."""
        return self._response_dict

    def _is_server_ready(self):
        """
        Check if the inference server is live and accessible.

        This method checks the availability of the inference server and raises an exception if it is not live or
        accessible. It ensures that the server is properly running and can be used for inference with the Koina
        model. Note: This method is primarily for internal use and typically called during model initialization.

        :raises ValueError: If the server responds with a not live status
        :raises InferenceServerException: If an exception occured while querying the server for its status.
        """
        try:
            if not self.client.is_server_live():
                raise ValueError("Server not yet started.")
        except InferenceServerException as e:
            if self.url == "koina.proteomicsdb.org:443":
                if self.ssl:
                    raise InferenceServerException(
                        "The public koina network seems to be inaccessible at the moment. "
                        "Please notify ludwig.lautenbacher@tum.de."
                    ) from e
                else:
                    raise InferenceServerException("To use the public koina network you need to set `ssl=True`.") from e
            raise InferenceServerException("Unknown error occured.", e.status(), e.debug_details()) from e

    def _is_model_ready(self):
        """
        Check if the specified model is available on the server.

        This method checks if the specified Koina model is available on the inference server. If the model is not
        available, it raises an exception indicating that the model is not accessible at the provided server URL.
        Note: This method is primarily for internal use and typically called during model initialization.

        :raises ValueError: If the specified model is not available at the server.
        :raises InferenceServerException: If an exception occured while querying the server for available models.
        """
        try:
            if not self.client.is_model_ready(self.model_name):
                raise ValueError(f"The model {self.model_name} is not available at {self.url}")
        except InferenceServerException as e:
            raise InferenceServerException("Unknown error occured.", e.status(), e.debug_details()) from e

    def __get_inputs(self):
        """
        Retrieve the input names and datatypes for the model.

        This method fetches the names and data types of the input tensors for the Koina model and stores them in
        the 'model_inputs' attribute. Note: This method is for internal use and is typically called during model
        initialization.

        :raises InferenceServerException: If an exception occured while querying the server for model inputs.
        """
        try:
            self.model_inputs = {i.name: i.datatype for i in self.client.get_model_metadata(self.model_name).inputs}
        except InferenceServerException as e:
            raise InferenceServerException("Unknown error occured.", e.status(), e.debug_details()) from e

    def __get_outputs(self, targets: Optional[List] = None):
        """
        Retrieve the output names and datatypes for the model.

        This method fetches the names and data types of the output tensors for the Koina model and stores them in
        the 'model_outputs' attribute. If a list of target names is supplied, the tensors are filtered for those.
        In case that the targets contain a name that is not a valid output of the requested model, a ValueError is
        raised. Note: This method is for internal use and is typically called during model initialization.

        :param targets: An optional list of target names to filter the predictions for. If this is None, all targets
            are added to list of output tensors to predict.
        :raises ValueError: If a target supplied is not a valid output name of the requested model.
        :raises InferenceServerException: If an exception occured while querying the server for model metadata.

        """
        try:
            model_outputs = self.client.get_model_metadata(self.model_name).outputs
            model_targets = [out.name for out in model_outputs]

            if targets is None:
                targets = model_targets
            else:
                for target in targets:
                    if target not in model_targets:
                        raise ValueError(
                            f"The supplied target {target} is not a valid output target of the model. "
                            f"Valid targets are {model_targets}."
                        )
            for i in model_outputs:
                if i.name in targets:
                    self.model_outputs[i.name] = i.datatype
        except InferenceServerException as e:
            raise InferenceServerException("Unknown error occured.", e.status(), e.debug_details()) from e

    def __get_batchsize(self):
        """
        Get the maximum batch size supported by the model's configuration.

        This method determines the maximum batch size supported by the Koina model's configuration and stores it
        in the 'batchsize' attribute. Note: This method is for internal use and is typically called during model
        initialization.
        :raises InferenceServerException: If an exception occured while querying the server for the max batchsize.
        """
        try:
            self.batchsize = self.client.get_model_config(self.model_name).config.max_batch_size
        except InferenceServerException as e:
            raise InferenceServerException("Unknown error occured.", e.status(), e.debug_details()) from e

    @staticmethod
    def __get_batch_outputs(names: KeysView[str]) -> List[InferRequestedOutput]:
        """
        Create InferRequestedOutput objects for the given output names.

        This method generates InferRequestedOutput objects for the specified output names. InferRequestedOutput objects
        are used to request specific outputs when performing inference. Note: This method is for internal use and is
        typically called during inference.

        :param names: A list of output names for which InferRequestedOutput objects should be created.

        :return: A list of InferRequestedOutput objects.
        """
        return [InferRequestedOutput(name) for name in names]

    def __get_batch_inputs(self, data: Dict[str, np.ndarray]) -> List[InferInput]:
        """
        Prepare a list of InferInput objects for the input data.

        This method prepares a list of InferInput objects for the provided input data. InferInput objects are used to
        specify the input tensors and their data when performing inference. Note: This method is for internal use and
        is typically called during inference.

        :param data: A dictionary containing input data for inference. Keys are input names, and values are numpy arrays.

        :return: A list of InferInput objects for the input data.
        """
        batch_inputs = []
        for iname, idtype in self.model_inputs.items():
            batch_inputs.append(InferInput(iname, (len(data[next(iter(data))]), 1), idtype))
            batch_inputs[-1].set_data_from_numpy(data[iname].reshape(-1, 1).astype(self.type_convert[idtype]))
        return batch_inputs

    def __extract_predictions(self, infer_result: InferResult) -> Dict[str, np.ndarray]:
        """
        Extract the predictions from an inference result.

        This method extracts the predictions from an inference result and organizes them in a dictionary with output
        names as keys and corresponding arrays as values. Note: This method is for internal use and is typically called
        during inference.

        :param infer_result: The result of an inference operation.

        :return: A dictionary containing the extracted predictions. Keys are output names, and values are numpy arrays.
        """
        predictions = {}
        for oname in self.model_outputs.keys():
            predictions[oname] = infer_result.as_numpy(oname)
        return predictions

    def __predict_batch(self, data: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
        """
        Perform batch inference and return the predictions.

        This method performs batch inference on the provided input data using the configured Koina model and returns the
        predictions. Note: This method is for internal use and is typically called during inference.

        :param data: A dictionary containing input data for batch inference. Keys are input names, and values are numpy arrays.

        :return: A dictionary containing the model's predictions. Keys are output names, and values are numpy arrays
            representing the model's output.
        """
        batch_outputs = self.__get_batch_outputs(self.model_outputs.keys())
        batch_inputs = self.__get_batch_inputs(data)
        infer_result = self.client.infer(self.model_name, inputs=batch_inputs, outputs=batch_outputs)

        return self.__extract_predictions(infer_result)

    def __predict_sequential(self, data: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
        """
        Perform sequential inference and return the predictions.

        This method performs sequential inference on the provided input data using the configured Koina model. It processes
        the input data batch by batch and returns the predictions.
        Note: This method is for internal use and is typically called during inference.

        :param data: A dictionary containing input data for inference. Keys are input names, and values are numpy arrays.

        :return: A dictionary containing the model's predictions. Keys are output names, and values are numpy arrays representing
            the model's output.
        """
        predictions: Dict[str, np.ndarray] = {}
        for data_batch in tqdm(self.__slice_dict(data, self.batchsize), desc="Getting predictions"):
            pred_batch = self.__predict_batch(data_batch)
            if predictions:
                predictions = self.__merge_array_dict(predictions, pred_batch)
            else:
                predictions = pred_batch  # Only first iteration to initialize dict keys
        return predictions

    @staticmethod
    def __slice_dict(data: Dict[str, np.ndarray], batchsize: int) -> Generator[Dict[str, np.ndarray], None, None]:
        """
        Slice the input data into batches of a specified batch size.

        This method takes the input data and divides it into smaller batches, each containing 'batchsize' elements. It yields
        these batches one at a time, allowing for batched processing of input data. Note: This method is for internal use and
        is typically called during batched inference.

        :param data: A dictionary containing input data for batch inference. Keys are input names, and values are numpy arrays.
        :param batchsize: The desired batch size for slicing the data.

        :yield: A dictionary containing a batch of input data with keys and values corresponding to the input names and
            batched arrays.
        """
        len_inputs = list(data.values())[0].shape[0]
        for i in range(0, len_inputs, batchsize):
            dict_slice = {}
            for k, v in data.items():
                dict_slice[k] = v[i : i + batchsize]
            yield dict_slice

    @staticmethod
    def __merge_array_dict(d1: Dict[str, np.ndarray], d2: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
        """
        Merge two dictionaries of arrays.

        This method takes two dictionaries, 'd1' and 'd2', each containing arrays with identical keys. It merges the
        arrays from both dictionaries, creating a new dictionary with the same keys and combined arrays. Note: This
        method is for internal use and is typically called during batched inference.

        :param d1: A dictionary containing arrays.
        :param d2: Another dictionary containing arrays with the same keys as d1.

        :raises NotImplementedError: If the keys in 'd1' and 'd2' do not match.
        :return: A dictionary containing merged arrays with the same keys as d1 and d2.

        Example:
        ```
        dict1 = {"output1": np.array([1.0, 2.0, 3.0]), "output2": np.array([4.0, 5.0, 6.0])}
        dict2 = {"output1": np.array([7.0, 8.0, 9.0]), "output2": np.array([10.0, 11.0, 12.0])}
        merged_dict = model.__merge_array_dict(dict1, dict2)
        print(merged_dict)
        ```
        """
        if d1.keys() != d2.keys():
            raise NotImplementedError(f"Keys in dictionary need to be equal {d1.keys(), d2.keys()}")
        out = {}
        for k in d1.keys():
            out[k] = np.concatenate([d1[k], d2[k]])
        return out

    @staticmethod
    def __merge_list_dict_array(dict_list: List[Dict[str, np.ndarray]]) -> Dict[str, np.ndarray]:
        """
        Merge a list of dictionaries of arrays.

        This method takes a list of dictionaries, where each dictionary contains arrays with identical keys. It merges
        the arrays from all dictionaries in the list, creating a new dictionary with the same keys and combined arrays.
        Note: This method is for internal use and is typically called during batched inference.

        :param dict_list: A list of dictionaries, each containing arrays with the same keys.
        :raises NotImplementedError: If the keys of all dictionaries in the list do not match.

        :return: A dictionary containing merged arrays with the same keys as the dictionaries in the list.

        Example::
            dict_list = [
                {"output1": np.array([1.0, 2.0, 3.0]), "output2": np.array([4.0, 5.0, 6.0])},
                {"output1": np.array([7.0, 8.0, 9.0]), "output2": np.array([10.0, 11.0, 12.0])},
                {"output1": np.array([13.0, 14.0, 15.0]), "output2": np.array([16.0, 17.0, 18.0])},
            ]
            merged_dict = model.__merge_list_dict_array(dict_list)
            print(merged_dict)
        """
        tmp = [x.keys() for x in dict_list]
        if not np.all([tmp[0] == x for x in tmp]):
            raise NotImplementedError(f"Keys of all dictionaries in the list need to be equal {tmp}")
        out = {}
        for k in tmp[0]:
            out[k] = np.concatenate([x[k] for x in dict_list])
        return out

    def __async_callback(
        self,
        infer_results: Dict[int, Union[Dict[str, np.ndarray], InferenceServerException]],
        request_id: int,
        result: Optional[InferResult],
        error: Optional[InferenceServerException],
    ):
        """
        Callback function for asynchronous inference.

        This method serves as a callback function for asynchronous inference. It is invoked when an asynchronous
        inference task is completed. The result of the task is appended to the 'infer_results' list, and any
        encountered error is checked and handled appropriately. Note: This method is for internal use and is typically
        called during asynchronous inference.

        :param infer_results: A dictionary to which the results of asynchronous inference will be added.
        :param request_id: The request id used as key in the infer_results dictionary
        :param result: The result of an asynchronous inference operation.
        :param error: An error, if any, encountered during asynchronous inference.
        """
        if error:
            infer_results[request_id] = error
        else:
            infer_results[request_id] = self.__extract_predictions(result)

    def __async_predict_batch(
        self,
        data: Dict[str, np.ndarray],
        infer_results: Dict[int, Union[Dict[str, np.ndarray], InferenceServerException]],
        request_id: int,
        timeout: int = 10000,
        retries: int = 10,
    ):
        """
        Perform asynchronous batch inference on the given data using the Koina model.

        This method initiates asynchronous batch inference on the provided input data using the configured Koina model.
        Results will be appended to the 'infer_results' list as they become available. The 'id' parameter is used to
        identify and order the results. The method will return when the inference request is completed or when the
        'timeout' is reached.

        :param data: A dictionary containing input data for batch inference. Keys are input names, and values are numpy arrays.
        :param infer_results: A dictionary to which the results of asynchronous inference will be added.
        :param request_id: An identifier for the inference request, used to track the order of completion.
        :param timeout: The maximum time (in seconds) to wait for the inference to complete. Defaults to 10 seconds.
        :param retries: The maximum number of requests in case of failure
        :yield: None, this is to separate async clien infer from checking the result
        """
        batch_outputs = self.__get_batch_outputs(self.model_outputs.keys())
        batch_inputs = self.__get_batch_inputs(data)

        for i in range(retries):
            if i > 0:  # need to yield first, before doing sth, but only after first time
                yield
                if isinstance(infer_results.get(request_id), InferResult):
                    break
                del infer_results[request_id]  # avoid race condition in case inference is slower than tqdm loop

            self.client.async_infer(
                model_name=self.model_name,
                request_id=str(request_id),
                inputs=batch_inputs,
                callback=partial(self.__async_callback, infer_results, request_id),
                outputs=batch_outputs,
                client_timeout=timeout,
            )

    def predict(
        self,
        data: Union[Dict[str, np.ndarray], pd.DataFrame],
        _async: bool = True,
        debug=False,
    ) -> Dict[str, np.ndarray]:
        """
        Perform inference on the given data using the Koina model.

        This method allows you to perform inference on the provided input data using the configured Koina model. You can
        choose to perform inference asynchronously (in parallel) or sequentially, depending on the value of the '_async'
        parameter. If asynchronous inference is selected, the method will return when all inference tasks are complete.
        Note: Ensure that the model and server are properly configured and that the input data matches the model's
        nput requirements.

        :param data: A dictionary or dataframe containing input data for inference. For the dictionary, keys are input names,
            and values are numpy arrays. In case of a dataframe, the input fields for the requested model must be present
            in the column names.
        :param _async: If True, perform asynchronous inference; if False, perform sequential inference. Defaults to True.
        :param debug: If True and using _async mode, store raw InferResult / InferServerException dictionary for later analysis.

        :return: A dictionary containing the model's predictions. Keys are output names, and values are numpy arrays
            representing the model's output.

        Example::
            model = Koina("Prosit_2019_intensity")
            input_data = {
                "peptide_sequences": np.array(["PEPTIDEK" for _ in range(size)]),
                "precursor_charges": np.array([2 for _ in range(size)]),
                "collision_energies": np.array([20 for _ in range(size)]),
                "fragmentation_types": np.array(["HCD" for _ in range(size)]),
                "instrument_types": np.array(["QE" for _ in range(size)])
            }
            predictions = model.predict(input_data)
        """
        if isinstance(data, pd.DataFrame):
            data = {
                input_field: data[alternative_column_map[input_field]].to_numpy()
                for input_field in self.model_inputs.keys()
            }
        if _async:
            return self.__predict_async(data, debug=debug)
        else:
            return self.__predict_sequential(data)

    def __predict_async(self, data: Dict[str, np.ndarray], debug=False) -> Dict[str, np.ndarray]:
        """
        Perform asynchronous inference on the given data using the Koina model.

        This method performs asynchronous inference on the provided input data using the configured Koina model.
        Asynchronous inference allows for parallel processing of input data, potentially leading to faster results.
        The method will return when all asynchronous inference tasks are complete. Note: Ensure that the model and server
        are properly configured and that the input data matches the model's input requirements.

        :param data: A dictionary containing input data for inference. Keys are input names, and values are numpy arrays.
        :param debug: If True, store raw InferResult / InferServerException dictionary for later analysis.

        :return: A dictionary containing the model's predictions. Keys are output names, and values are numpy arrays
            representing the model's output.
        """
        infer_results: Dict[int, Union[Dict[str, np.ndarray], InferenceServerException]] = {}
        tasks = []
        for i, data_batch in enumerate(self.__slice_dict(data, self.batchsize)):
            tasks.append(self.__async_predict_batch(data_batch, infer_results, request_id=i, retries=3))
            next(tasks[i])

        n_tasks = i + 1
        with tqdm(total=n_tasks, desc="Getting predictions", disable=self.disable_progress_bar) as pbar:
            unfinished_tasks = [i for i in range(n_tasks)]
            while pbar.n < n_tasks:
                time.sleep(0.5)
                new_unfinished_tasks = []
                for j in unfinished_tasks:
                    result = infer_results.get(j)
                    if result is None:
                        new_unfinished_tasks.append(j)
                    elif isinstance(result, dict):
                        pbar.n += 1
                    else:  # unexpected result / exception -> try again
                        try:
                            next(tasks[j])
                            new_unfinished_tasks.append(j)
                        except StopIteration:
                            pbar.n += 1

                unfinished_tasks = new_unfinished_tasks
                pbar.refresh()

        return self.__handle_results(infer_results, debug)

    def __handle_results(
        self, infer_results: Dict[int, Union[Dict[str, np.ndarray], InferenceServerException]], debug: bool
    ) -> Dict[str, np.ndarray]:
        """
        Handles the results.

        :param infer_results: The dictionary containing the inferred results
        :param debug: whether to store the infer_results in the response_dict attribute

        :raises InferenceServerException: If at least one batch of predictions could not be inferred.

        :return: A dictionary containing the model's predictions. Keys are output names, and values are numpy arrays
            representing the model's output.
        """
        if debug:
            self._response_dict = infer_results
        try:
            # sort according to request id
            infer_results_to_return = [infer_results[i] for i in range(len(infer_results))]
            return self.__merge_list_dict_array(infer_results_to_return)
        except AttributeError:
            for res in infer_results.values():
                if isinstance(res, InferenceServerException):
                    warnings.warn(res.message(), stacklevel=1)
            else:
                raise InferenceServerException(
                    """
                    At least one request failed. Check the error message above and try again.
                    To get a list of responses run koina.predict(..., debug = True), then call koina.response_dict
                    """
                ) from None
