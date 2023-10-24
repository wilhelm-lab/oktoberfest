import logging
import re
from math import ceil
from multiprocessing import current_process
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from spectrum_fundamentals.metrics.similarity import SimilarityMetrics
from tqdm.auto import tqdm
from tritonclient.grpc import InferenceServerClient, InferInput, InferRequestedOutput

from ..data.spectra import FragmentType, Spectra

logger = logging.getLogger(__name__)


def grpc_predict(
    library: Spectra,
    url: str,
    intensity_model: str,
    irt_model: str,
    ssl: bool = True,
    alignment: bool = False,
    job_type: str = "",
):
    """
    Use grpc to predict library and add predictions to library.

    :param library: Spectra object with the library
    :param url: Url including the port of the prediction server
    :param intensity_model: the name of the intensity model on the server
    :param irt_model: the name of the irt model on the server
    :param ssl: whether or not the server requires an ssl encrypted transportation, default = True
    :param alignment: True if alignment present
    :param job_type: TODO
    :return: grpc predictions if we are trying to generate spectral library
    """
    triton_client = InferenceServerClient(url=url, ssl=ssl)
    batch_size = 1000

    intensity_outputs = ["intensities", "mz", "annotation"]
    intensity_input_data = {
        "peptide_sequences": (
            library.spectra_data["MODIFIED_SEQUENCE"].to_numpy().reshape(-1, 1).astype(np.object_),
            "BYTES",
        ),
        "collision_energies": (
            library.spectra_data["COLLISION_ENERGY"].to_numpy().reshape(-1, 1).astype(np.float32),
            "FP32",
        ),
        "precursor_charges": (
            library.spectra_data["PRECURSOR_CHARGE"].to_numpy().reshape(-1, 1).astype(np.int32),
            "INT32",
        ),
    }
    if "tmt" in intensity_model.lower() or "ptm" in intensity_model.lower():
        intensity_input_data["fragmentation_types"] = (
            library.spectra_data["FRAGMENTATION"].to_numpy().reshape(-1, 1).astype(np.object_),
            "BYTES",
        )

    intensity_predictions = infer_predictions(
        triton_client,
        model=intensity_model,
        input_data=intensity_input_data,
        outputs=intensity_outputs,
        batch_size=batch_size,
    )
    intensity_predictions["intensities"][np.where(intensity_predictions["intensities"] < 1e-7)] = 0.0

    irt_input_data = {"peptide_sequences": intensity_input_data["peptide_sequences"]}
    irt_outputs = ["irt"]
    irt_predictions = infer_predictions(
        triton_client,
        model=irt_model,
        input_data=irt_input_data,
        outputs=irt_outputs,
        batch_size=batch_size,
    )

    if job_type == "SpectralLibraryGeneration":
        intensity_prediction_dict = {
            "intensity": intensity_predictions["intensities"],
            "fragmentmz": intensity_predictions["mz"],
            "annotation": parse_fragment_labels(
                intensity_predictions["annotation"],
                library.spectra_data["PRECURSOR_CHARGE"].to_numpy()[:, None],
                library.spectra_data["PEPTIDE_LENGTH"].to_numpy()[:, None],
            ),
        }
        output_dict = {intensity_model: intensity_prediction_dict, irt_model: irt_predictions["irt"]}
        return output_dict

    intensities_pred = pd.DataFrame()
    intensities_pred["intensity"] = intensity_predictions["intensities"].tolist()
    library.add_matrix(intensities_pred["intensity"], FragmentType.PRED)

    if alignment:
        return

    library.add_column(irt_predictions["irt"], name="PREDICTED_IRT")


def infer_predictions(
    triton_client: InferenceServerClient,
    model: str,
    input_data: Dict[str, Tuple[np.ndarray, str]],
    outputs: List[str],
    batch_size: int,
) -> Dict[str, np.ndarray]:
    """
    Infer predictions from a triton client.

    :param triton_client: An inference client using grpc
    :param model: a model that is recognized by the server specified in the triton client
    :param input_data: a dictionary that contains the input names (key) for the specific model
        and a tuple of the input_data as a numpy array of shape [:, 1] and the dtype recognized
        by the triton client (value).
    :param outputs: a list of output names for the specific model
    :param batch_size: the number of elements from the input_data that should be provided to the
        triton client at once
    :return: a dictionary containing the predictions (values) for the given outputs (keys)
    """
    num_spec = len(input_data[list(input_data)[0]][0])
    predictions: Dict[str, List[np.ndarray]] = {output: [] for output in outputs}

    n_batches = ceil(num_spec / batch_size)
    process_identity = current_process()._identity
    if len(process_identity) > 0:
        position = process_identity[0]
    else:
        position = 0

    with tqdm(
        total=n_batches,
        position=position,
        desc=f"Inferring predictions for {num_spec} spectra with batch site {batch_size}",
        leave=True,
    ) as progress:
        for i in range(0, n_batches):
            progress.update(1)
            # logger.info(f"Predicting batch {i+1}/{n_batches}.")
            infer_inputs = []
            for input_key, (data, dtype) in input_data.items():
                batch_data = data[i * batch_size : (i + 1) * batch_size]
                infer_input = InferInput(input_key, batch_data.shape, dtype)
                infer_input.set_data_from_numpy(batch_data)
                infer_inputs.append(infer_input)

            infer_outputs = [InferRequestedOutput(output) for output in outputs]

            prediction = triton_client.infer(model, inputs=infer_inputs, outputs=infer_outputs)

            for output in outputs:
                predictions[output].append(prediction.as_numpy(output))

    return {key: np.vstack(value) for key, value in predictions.items()}


def parse_fragment_labels(spectra_labels: np.ndarray, precursor_charges: np.ndarray, seq_lengths: np.ndarray):
    """Uses regex to parse labels."""
    pattern = rb"([y|b])([0-9]{1,2})\+([1-3])"
    fragment_types = []
    fragment_numbers = []
    fragment_charges = []
    for spectrum_labels in spectra_labels:
        types = []
        numbers = []
        charges = []
        for label in spectrum_labels:
            match = re.match(pattern, label)
            if match:
                groups = match.groups()
                types.append(groups[0].decode())
                numbers.append(int(groups[1]))
                charges.append(int(groups[2]))
            else:
                raise ValueError(f"String {label} does not match the expected fragment label pattern")
        fragment_types.append(types)
        fragment_numbers.append(numbers)
        fragment_charges.append(charges)

    fragment_type_array = np.array(fragment_types)
    fragment_number_array = np.array(fragment_numbers)
    fragment_charge_array = np.array(fragment_charges)
    mask = np.where((fragment_charge_array > precursor_charges) | (fragment_number_array >= seq_lengths))
    fragment_type_array[mask] = "N"
    fragment_number_array[mask] = 0
    fragment_charge_array[mask] = 0

    return {"type": fragment_type_array, "number": fragment_number_array, "charge": fragment_charge_array}


def _prepare_alignment_df(library: Spectra) -> Spectra:
    """
    Prepare an alignment DataFrame from the given Spectra library.

    This function creates an alignment DataFrame by removing decoy and HCD fragmented spectra
    from the input library, selecting the top 1000 highest-scoring spectra, and repeating the
    DataFrame for each collision energy (CE) in the range [18, 50].

    :param library: the library to be propagated
    :return: a library that is modified according to the description above
    """
    alignment_library = Spectra()
    alignment_library.spectra_data = library.spectra_data.copy()

    # Remove decoy and HCD fragmented spectra
    alignment_library.spectra_data = alignment_library.spectra_data[
        (alignment_library.spectra_data["FRAGMENTATION"] == "HCD") & (~alignment_library.spectra_data["REVERSE"])
    ]
    # Select the 1000 highest scoring or all if there are less than 1000
    alignment_library.spectra_data = alignment_library.spectra_data.sort_values(by="SCORE", ascending=False).iloc[:1000]

    # Repeat dataframe for each CE
    ce_range = range(18, 50)
    nrow = len(alignment_library.spectra_data)
    alignment_library.spectra_data = pd.concat([alignment_library.spectra_data for _ in ce_range], axis=0)
    alignment_library.spectra_data["COLLISION_ENERGY"] = np.repeat(ce_range, nrow)
    alignment_library.spectra_data.reset_index(inplace=True)
    return alignment_library


def ce_calibration(library: Spectra, **server_kwargs) -> pd.Series:
    """
    Calculate best collision energy for peptide property predictions.

    The function propagates the provided library object to test NCEs in th range [18,50], performs
    intensity prediction for the 1000 highest scoring target PSMs at each NCE and computes the spectral angle
    between predicted and observed intensities before returning the alignment library.

    :param library: spectral library to perform CE calibration on
    :param server_kwargs: Additional parameters that are forwarded to grpc_predict
    :return: pandas series containing the spectral angle for all tested collision energies
    """
    alignment_library = _prepare_alignment_df(library)
    grpc_predict(alignment_library, alignment=True, **server_kwargs)
    _alignment(alignment_library)
    return alignment_library


def _alignment(alignment_library: Spectra):
    """
    Perform the alignment of predicted versus raw intensities.

    The function calculates the spectral angle between predicted and observed fragment intensities and
    adds it as a column to the alignment library.

    :param alignment_library: the library to perform the alignment on
    """
    pred_intensity = alignment_library.get_matrix(FragmentType.PRED)
    raw_intensity = alignment_library.get_matrix(FragmentType.RAW)
    # return pred_intensity.toarray(), raw_intensity.toarray()
    sm = SimilarityMetrics(pred_intensity, raw_intensity)
    alignment_library.spectra_data["SPECTRAL_ANGLE"] = sm.spectral_angle(raw_intensity, pred_intensity, 0)

class KoinaModel:
    """
    A class for interacting with Koina models for inference.

    Parameters:
    model_name (str): The name of the Koina model.
    server_url (str, optional): The URL of the Inference Server. Defaults to "koina.proteomicsdb.org:443".
    ssl (bool, optional): Use SSL for communication with the server. Defaults to True.

    Methods:
    - predict(data, disable_progress_bar, async_): Perform inference (sequential or asynchronous).
    - _is_server_ready(): Check if the Inference Server is live and accessible.
    - _is_model_ready(): Check if the specified model is available on the server.
    - __get_inputs(): Retrieve the input names and datatypes for the model.
    - __get_outputs(): Retrieve the output names and datatypes for the model.
    - __get_batchsize(): Get the maximum batch size supported by the model's configuration.
    - __get_batch_outputs(names): Create InferRequestedOutput objects for the given output names.
    - __get_batch_inputs(data): Prepare a list of InferInput objects for the input data.
    - __extract_predictions(infer_result): Extract the predictions from an inference result.
    - __predict_batch(data): Perform batch inference and return the predictions.
    - __predict_sequential(data, disable_progress_bar): Perform sequential inference and return the predictions.
    - __slice_dict(data, batchsize): Slice the input data into batches of a specified batch size.
    - __merge_array_dict(d1, d2): Merge two dictionaries of arrays.
    - __merge_list_dict_array(list): Merge a list of dictionaries of arrays.
    - __async_callback(infer_results, result, error): Callback function for asynchronous inference.
    - __async_predict_batch(data, infer_results, id, timeout): Perform asynchronous batch inference.
    - __predict_async(data, disable_progress_bar): Perform asynchronous inference and return the predictions.
    """

    def __init__(self, model_name, server_url="koina.proteomicsdb.org:443", ssl=True):
        """
        Initialize a KoinaModel instance with the specified parameters.

        Parameters:
        - model_name (str): The name of the Koina model to be used for inference.
        - server_url (str, optional): The URL of the Inference Server. Defaults to "koina.proteomicsdb.org:443".
        - ssl (bool, optional): Indicates whether to use SSL for communication with the server. Defaults to True.

        This constructor initializes the KoinaModel instance, connecting it to the specified Inference Server. It checks the availability of the server, the specified model, retrieves input and output information, and determines the maximum batch size supported by the model's configuration.

        Note: To use this class, ensure that the Inference Server is properly configured and running, and that the specified model is available on the server.
        """
        self.model_inputs = {}
        self.model_outputs = {}
        self.batchsize = None
        
        self.model_name = model_name
        self.url = server_url
        self.ssl = ssl
        self.client = InferenceServerClient(
            url=server_url, ssl=ssl
        )

        self.type_convert = {
            'FP32': np.dtype('float32'),
            'BYTES': np.dtype('O'),
            'INT16': np.dtype('int16'),
            'INT32': np.dtype('int32'),
            'INT64':np.dtype('int64'),
        }

        self._is_server_ready()
        self._is_model_ready()

        self.__get_inputs()
        self.__get_outputs()
        self.__get_batchsize()

    
    def _is_server_ready(self):
        """
        Check if the Inference Server is live and accessible.

        This method checks the availability of the Inference Server and raises an exception if it is not live or accessible. It ensures that the server is properly running and can be used for inference with the Koina model.

        Note: This method is primarily for internal use and typically called during model initialization.

        Raises:
        - ValueError: If the server is not live or accessible.
        """
        try:
            if not self.client.is_server_live():
                raise ValueError("Server not yet startet")
        except InferenceServerException as e:
            if self.url == "koina.proteomicsdb.org:443":
                if self.ssl:
                    raise ValueError("The public koina network seems to be inaccessible at the moment please notify Ludwig.Lautenbacher@tum.de.")
                else:
                    raise ValueError("The use the public koina network you need to set `ssl=True`.")
            else:
                raise e

    def _is_model_ready(self):
        """
        Check if the specified model is available on the server.

        This method checks if the specified Koina model is available on the Inference Server. If the model is not available, it raises an exception indicating that the model is not accessible at the provided server URL.

        Note: This method is primarily for internal use and typically called during model initialization.

        Raises:
        - ValueError: If the specified model is not available at the server.
        """
        if not self.client.is_model_ready(self.model_name):
            raise ValueError(f"The model {self.model_name} is not available at {self.url}")
        
    def __get_inputs(self):
        """
        Retrieve the input names and datatypes for the model.

        This method fetches the names and data types of the input tensors for the Koina model and stores them in the 'model_inputs' attribute.

        Note: This method is for internal use and is typically called during model initialization.
        """
        for i in self.client.get_model_metadata(self.model_name).inputs:
            self.model_inputs[i.name] = i.datatype

    def __get_outputs(self):
        """
        Retrieve the output names and datatypes for the model.

        This method fetches the names and data types of the output tensors for the Koina model and stores them in the 'model_outputs' attribute.

        Note: This method is for internal use and is typically called during model initialization.
        """
        for i in self.client.get_model_metadata(self.model_name).outputs:
            self.model_outputs[i.name] = i.datatype

    def __get_batchsize(self):
        """
        Get the maximum batch size supported by the model's configuration.

        This method determines the maximum batch size supported by the Koina model's configuration and stores it in the 'batchsize' attribute.

        Note: This method is for internal use and is typically called during model initialization.
        """
        self.batchsize = self.client.get_model_config(self.model_name).config.max_batch_size

    @staticmethod
    def __get_batch_outputs(names):
        """
        Create InferRequestedOutput objects for the given output names.

        Parameters:
        - names (list): A list of output names for which InferRequestedOutput objects should be created.

        Returns:
        - list: A list of InferRequestedOutput objects.

        This method generates InferRequestedOutput objects for the specified output names. InferRequestedOutput objects are used to request specific outputs when performing inference.

        Note: This method is for internal use and is typically called during inference.
        """
        return [InferRequestedOutput(name) for name in names]

    def __get_batch_inputs(self,data):
        """
        Prepare a list of InferInput objects for the input data.

        Parameters:
        - data (dict): A dictionary containing input data for inference. Keys are input names, and values are numpy arrays.

        Returns:
        - list: A list of InferInput objects for the input data.

        This method prepares a list of InferInput objects for the provided input data. InferInput objects are used to specify the input tensors and their data when performing inference.

        Note: This method is for internal use and is typically called during inference.
        """
        batch_inputs = []
        for iname, idtype in self.model_inputs.items():
            batch_inputs.append(
                InferInput(iname, (len(data[next(iter(data))]),1), idtype)
            )
            batch_inputs[-1].set_data_from_numpy(data[iname].reshape(-1,1).astype(self.type_convert[idtype]))
        return batch_inputs

    def __extract_predictions(self, infer_result):
        """
        Extract the predictions from an inference result.

        Parameters:
        - infer_result: The result of an inference operation.

        Returns:
        - dict: A dictionary containing the extracted predictions. Keys are output names, and values are numpy arrays.

        This method extracts the predictions from an inference result and organizes them in a dictionary with output names as keys and corresponding arrays as values.

        Note: This method is for internal use and is typically called during inference.
        """
        predictions = {}
        for oname in self.model_outputs.keys():
            predictions[oname] = infer_result.as_numpy(oname)
        return predictions

    def __predict_batch(self, data):
        """
        Perform batch inference and return the predictions.

        Parameters:
        - data (dict): A dictionary containing input data for batch inference. Keys are input names, and values are numpy arrays.

        Returns:
        - dict: A dictionary containing the model's predictions. Keys are output names, and values are numpy arrays representing the model's output.

        This method performs batch inference on the provided input data using the configured Koina model and returns the predictions.

        Note: This method is for internal use and is typically called during inference.
        """
        batch_outputs = self.__get_batch_outputs(self.model_outputs.keys())
        batch_inputs = self.__get_batch_inputs(data)

        infer_result =  self.client.infer(self.model_name, inputs=batch_inputs, outputs=batch_outputs)

        return self.__extract_predictions(infer_result)
    
    def __predict_sequential(self, data, disable_progress_bar):
        """
        Perform sequential inference and return the predictions.

        Parameters:
        - data (dict): A dictionary containing input data for inference. Keys are input names, and values are numpy arrays.
        - disable_progress_bar (bool): If True, disable the progress bar during inference.

        Returns:
        - dict: A dictionary containing the model's predictions. Keys are output names, and values are numpy arrays representing the model's output.

        This method performs sequential inference on the provided input data using the configured Koina model. It processes the input data batch by batch and returns the predictions. You can choose to disable the progress bar during inference using the 'disable_progress_bar' parameter.

        Note: This method is for internal use and is typically called during inference.
        """
        predictions = {}
        for data_batch in tqdm(self.__slice_dict(data, self.batchsize), desc="Getting predictions", disable=disable_progress_bar):
            pred_batch = self.__predict_batch(data_batch)
            if predictions:
                predictions = self.__merge_array_dict(predictions,pred_batch)
            else:
                predictions = pred_batch # Only first iteration to initialize dict keys
        return predictions

    @staticmethod
    def __slice_dict(data, batchsize):
        """
        Slice the input data into batches of a specified batch size.

        Parameters:
        - data (dict): A dictionary containing input data for batch inference. Keys are input names, and values are numpy arrays.
        - batchsize (int): The desired batch size for slicing the data.

        Yields:
        - dict: A dictionary containing a batch of input data with keys and values corresponding to the input names and batched arrays.

        This method takes the input data and divides it into smaller batches, each containing 'batchsize' elements. It yields these batches one at a time, allowing for batched processing of input data.

        Note: This method is for internal use and is typically called during batched inference.
        """
        len_inputs = list(data.values())[0].shape[0]
        for i in range(0, len_inputs, batchsize):
            dict_slice = {}
            for k,v in data.items():
                dict_slice[k] = v[i:i+batchsize]
            yield dict_slice
    
    @staticmethod
    def __merge_array_dict(d1, d2):
        """
        Merge two dictionaries of arrays.

        Parameters:
        - d1 (dict): A dictionary containing arrays.
        - d2 (dict): Another dictionary containing arrays with the same keys as d1.

        Returns:
        - dict: A dictionary containing merged arrays with the same keys as d1 and d2.

        This method takes two dictionaries, 'd1' and 'd2', each containing arrays with identical keys. It merges the arrays from both dictionaries, creating a new dictionary with the same keys and combined arrays.

        Note: This method is for internal use and is typically called during batched inference.

        Raises:
        - NotImplementedError: If the keys in 'd1' and 'd2' do not match.

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
            out[k] = np.concatenate([d1[k],d2[k]])
        return out
    
    @staticmethod
    def __merge_list_dict_array(list):
        """
        Merge a list of dictionaries of arrays.

        Parameters:
        - list (list): A list of dictionaries, each containing arrays with the same keys.

        Returns:
        - dict: A dictionary containing merged arrays with the same keys as the dictionaries in the list.

        This method takes a list of dictionaries, where each dictionary contains arrays with identical keys. It merges the arrays from all dictionaries in the list, creating a new dictionary with the same keys and combined arrays.

        Note: This method is for internal use and is typically called during batched inference.

        Raises:
        - NotImplementedError: If the keys of all dictionaries in the list do not match.

        Example:
        ```
        list_of_dicts = [
            {"output1": np.array([1.0, 2.0, 3.0]), "output2": np.array([4.0, 5.0, 6.0])},
            {"output1": np.array([7.0, 8.0, 9.0]), "output2": np.array([10.0, 11.0, 12.0])},
            {"output1": np.array([13.0, 14.0, 15.0]), "output2": np.array([16.0, 17.0, 18.0])},
        ]
        merged_dict = model.__merge_list_dict_array(list_of_dicts)
        print(merged_dict)
        ```
        """
        tmp = [x.keys() for x in list]
        if not np.all([tmp[0] == x for x in tmp]):
            raise NotImplementedError(f"Keys of all dictionaries in the list need to be equal {tmp}")
        out = {}
        for k in tmp[0]:
            out[k] = np.concatenate([x[k] for x in list])
        return out

    def __async_callback(self, infer_results, result, error):
        """
        Callback function for asynchronous inference.

        Parameters:
        - infer_results (list): A list to which the results of asynchronous inference will be appended.
        - result: The result of an asynchronous inference operation.
        - error: An error, if any, encountered during asynchronous inference.

        This method serves as a callback function for asynchronous inference. It is invoked when an asynchronous inference task is completed. The result of the task is appended to the 'infer_results' list, and any encountered error is checked and handled appropriately.

        Note: This method is for internal use and is typically called during asynchronous inference.
        """
        if error:
            raise error
        else:
            infer_results.append(result)

    def __async_predict_batch(self, data, infer_results, id, timeout=10):
        """
        Perform asynchronous batch inference on the given data using the Koina model.

        Parameters:
        - data (dict): A dictionary containing input data for batch inference. Keys are input names, and values are numpy arrays.
        - infer_results (list): A list to which the results of asynchronous inference will be appended.
        - id (int): An identifier for the inference request, used to track the order of completion.
        - timeout (int, optional): The maximum time (in seconds) to wait for the inference to complete. Defaults to 10 seconds.

        This method initiates asynchronous batch inference on the provided input data using the configured Koina model. Results will be appended to the 'infer_results' list as they become available. The 'id' parameter is used to identify and order the results. The method will return when the inference request is completed or when the 'timeout' is reached.
        """
        batch_outputs = self.__get_batch_outputs(self.model_outputs.keys())
        batch_inputs = self.__get_batch_inputs(data)

        self.client.async_infer(
            model_name=self.model_name,
            request_id = str(id),
            inputs=batch_inputs,
            callback=partial(self.__async_callback, infer_results),
            outputs=batch_outputs,
            client_timeout=timeout,
        )

    def predict(self, data, disable_progress_bar=False, async_=True):
        """
        Perform inference on the given data using the Koina model.

        Parameters:
        - data (dict): A dictionary containing input data for inference. Keys are input names, and values are numpy arrays.
        - disable_progress_bar (bool, optional): If True, disable the progress bar during inference. Defaults to False.
        - async_ (bool, optional): If True, perform asynchronous inference; if False, perform sequential inference. Defaults to True.

        Returns:
        - dict: A dictionary containing the model's predictions. Keys are output names, and values are numpy arrays representing the model's output.

        This method allows you to perform inference on the provided input data using the configured Koina model. You can choose to perform inference asynchronously (in parallel) or sequentially, depending on the value of the 'async_' parameter. If asynchronous inference is selected, the method will return when all inference tasks are complete.

        Note: Ensure that the model and server are properly configured and that the input data matches the model's input requirements.

        Example:
        ```
        model = KoinaModel("Prosit_2019_intensity")
        input_data = {
            "peptide_sequences": np.array(["PEPTIDEK" for _ in range(size)]),
            "precursor_charges": np.array([2 for _ in range(size)]),
            "collision_energies": np.array([20 for _ in range(size)]),
            "fragmentation_types": np.array(["HCD" for _ in range(size)]),
            "instrument_types": np.array(["QE" for _ in range(size)])
        }
        predictions = model.predict(input_data)
        ```
        """
        if async_:
            self.__predict_async(data, disable_progress_bar=disable_progress_bar)
        else:
            self.__predict_sequential(data, disable_progress_bar=disable_progress_bar)

    def __predict_async(self, data, disable_progress_bar=False):
        """
        Perform asynchronous inference on the given data using the Koina model.

        Parameters:
        - data (dict): A dictionary containing input data for inference. Keys are input names, and values are numpy arrays.
        - disable_progress_bar (bool, optional): If True, disable the progress bar during asynchronous inference. Defaults to False.

        Returns:
        - dict: A dictionary containing the model's predictions. Keys are output names, and values are numpy arrays representing the model's output.

        This method performs asynchronous inference on the provided input data using the configured Koina model. Asynchronous inference allows for parallel processing of input data, potentially leading to faster results. The method will return when all asynchronous inference tasks are complete.

        Note: Ensure that the model and server are properly configured and that the input data matches the model's input requirements.
        """
        infer_results = []
        for i,data_batch in enumerate(self.__slice_dict(data, self.batchsize)):
            self.__async_predict_batch(data_batch, infer_results,id=i)
        
        with tqdm(total=i+1, desc="Getting predictions", disable=disable_progress_bar) as pbar:
            while len(infer_results) != i+1:
                pbar.n = len(infer_results)
                pbar.refresh()
                time.sleep(1)
            pbar.n = len(infer_results)
            pbar.refresh()

        # sort according to request id
        infer_results = [self.__extract_predictions(infer_results[i]) for i in np.argsort(np.array([int(y.get_response("id")["id"]) for y in infer_results]))]
        
        return self.__merge_list_dict_array(infer_results)