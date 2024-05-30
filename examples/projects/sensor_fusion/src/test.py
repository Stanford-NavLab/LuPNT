import numpy as np
from numba import njit

from pylupnt.crater_detection.matching.utils import shift_nd
from pylupnt.crater_detection.matching import CraterDatabase
import pylupnt as pnt

db = CraterDatabase.from_file(
    pnt.utils.find_file("lunar_crater_database_robbins_2018.csv"),
    latlims=(0, 30),
    longlims=(0, 30),
    diamlims=(4, 40),
)
print("Database loaded")

# def main():
#     # Create a sample array and shift values
#     sample_array = np.array(
#         [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]], dtype=np.float64
#     )
#     shift_values = np.array([1, 2, 3], dtype=np.int64)

#     # Print original array
#     print("Original array:")
#     print(sample_array)

#     # Print shift values
#     print("Shift values:")
#     print(shift_values)

#     # Apply shift_nd function
#     shifted_array = shift_nd(sample_array, shift_values)

#     # Print shifted array
#     print("Shifted array:")
#     print(shifted_array)


# if __name__ == "__main__":
#     main()
