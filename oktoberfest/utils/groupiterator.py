from typing import Optional

import pandas as pd


def group_iterator(df, group_by_column: str, max_batch_size: Optional[int] = None):
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
