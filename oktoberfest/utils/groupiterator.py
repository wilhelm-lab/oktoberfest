from typing import Optional, Union

import pandas as pd


def group_iterator(df: pd.DataFrame, group_by_column: Union[str, list[str]], max_batch_size: Optional[int] = None):
    """
    Returns an index iterator producing chunks for each group of a given max size.

    This function groups a given dataframe by the specified column(s) and provides an iterator
    yielding the dataframe indices in chunks of up to a given number of elements in each
    group. If the remaining elements of a group or the total number of elements in a group are
    less than the given number, the chunk is smaller, accordingly. If no upper limit is provided,
    the indices for the entire group are returned at once irrespective of the group size.

    :param df: The dataframe to produce the iterator from
    :param group_by_column: The name(s) of the columns to group the dataframe by
    :param max_batch_size: Optional upper limit for the number of indices within each group that
        are yieled at once.

    :yields: An iterator yielding the dataframe index in batches containing only one group and up
        to the provided number of elements at once
    """
    # Group the dataframe by the specified column
    grouped = df.groupby(group_by_column)

    for _, group_df in grouped:
        # Get the indices of the current group
        indices = group_df.index.to_numpy()
        if max_batch_size is None:
            yield indices
            continue

        # Calculate the number of full batches
        num_full_batches = len(indices) // max_batch_size

        # Yield full batches
        for i in range(num_full_batches):
            yield indices[i * max_batch_size : (i + 1) * max_batch_size]

        # Yield the remaining elements as the last batch
        remainder = len(indices) % max_batch_size
        if remainder > 0:
            yield indices[num_full_batches * max_batch_size :]
