#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DIPHA, PERSEUS <=> Numpy Format Converter

A utility script to convert between DIPHAformat and NumPy arrays.

Supported conversions:
- DIPHA complex files (.complex) to NumPy arrays (.npy)
- PERSEUS text files (.txt) to NumPy arrays (.npy)
- NumPy arrays (.npy) to DIPHA complex files (.complex)
- DIPHA output/diagram files (.output/.diagram) to NumPy arrays
- DiaMorse output text files (.txt) to NumPy arrays

"""

import struct
import numpy as np
import argparse
import os
from pathlib import Path
from typing import Tuple, Optional
import sys

def setup_argument_parser() -> argparse.ArgumentParser:
    """Set up command line argument parser with detailed help."""
    parser = argparse.ArgumentParser(
        prog="dipha2npy",
        description="Convert between DIPHA format and NumPy arrays for topological data analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python dipha2npy.py input.complex output.npy    # DIPHA to NumPy
  python dipha2npy.py input.npy output.complex    # NumPy to DIPHA
  python dipha2npy.py results.output analysis.npy # Persistence diagram to NumPy
  python dipha2npy.py diamorse.txt features.npy   # DiaMorse to NumPy
        """
    )

    parser.add_argument(
        'from_file',
        help='Input file path (supports: .complex, .npy, .output, .diagram, .txt)'
    )
    parser.add_argument(
        'to_file',
        help='Output file path'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output with detailed statistics'
    )

    return parser

def convert_dipha_complex_to_numpy(input_file: str, output_file: str, verbose: bool = False) -> None:
    """
    Convert DIPHA complex format to NumPy array.

    DIPHA complex files contain multidimensional arrays representing
    topological complexes used in persistent homology computations.

    Args:
        input_file: Path to input .complex file
        output_file: Path to output .npy file
        verbose: Enable detailed output
    """
    try:
        with open(input_file, 'rb') as f:
            data = f.read()

        # Read DIPHA header: magic number, type, size, dimensions
        magic, data_type, size, dimensions = struct.unpack_from("qqqq", data, 0)

        if verbose:
            print(f"DIPHA Magic: {magic}, Type: {data_type}, Size: {size}, Dimensions: {dimensions}")

        # Read shape information
        offset = 8 * 4  # 8 bytes per long, 4 longs in header
        shape = struct.unpack_from("q" * dimensions, data, offset)

        # Read the actual data
        data_offset = offset + 8 * dimensions
        array_data = struct.unpack_from("d" * size, data, data_offset)

        # Reshape the data
        result_array = np.array(array_data).reshape(shape)

        print(f"Converted complex to array: shape={result_array.shape}")
        print(f"Data range: [{result_array.min():.6f}, {result_array.max():.6f}]")
        print(f"Mean value: {result_array.mean():.6f}")

        np.save(output_file, result_array)
        print(f"Saved to: {output_file}")

    except Exception as e:
        print(f"Error converting DIPHA complex: {e}")
        sys.exit(1)

def convert_numpy_to_dipha_complex(input_file: str, output_file: str, verbose: bool = False) -> None:
    """
    Convert NumPy array to DIPHA complex format.

    This creates a DIPHA-compatible binary file that can be used
    with DIPHA for persistent homology computation.

    Args:
        input_file: Path to input .npy file
        output_file: Path to output .complex file
        verbose: Enable detailed output
    """
    try:
        data_array = np.load(input_file)

        print(f"Loaded NumPy array: shape={data_array.shape}")
        print(f"Data range: [{data_array.min():.6f}, {data_array.max():.6f}]")
        print(f"Mean value: {data_array.mean():.6f}")

        size = np.prod(data_array.shape)
        dimensions = len(data_array.shape)

        with open(output_file, "wb") as f:
            # Write DIPHA header
            # Magic number: 8067171840, Type: 1 (complex), Size, Dimensions
            f.write(struct.pack("q" * 4, 8067171840, 1, size, dimensions))

            # Write shape
            f.write(struct.pack("q" * dimensions, *data_array.shape))

            # Write data with appropriate transposition for DIPHA format
            if dimensions == 3:
                # Transpose for 3D: (z,y,x) ordering expected by DIPHA
                transposed_data = data_array.transpose((2, 1, 0))
            elif dimensions == 2:
                # Transpose for 2D: (y,x) ordering expected by DIPHA
                transposed_data = data_array.transpose((1, 0))
            else:
                transposed_data = data_array

            f.write(struct.pack("d" * size, *transposed_data.flatten()))

        print(f"Saved DIPHA complex to: {output_file}")

    except Exception as e:
        print(f"Error converting NumPy to DIPHA: {e}")
        sys.exit(1)

def convert_perseus_to_numpy(input_file: str, output_file: str, verbose: bool = False) -> None:
    """
    Convert PERSEUS text format to NumPy array.

    PERSEUS format structure:
    - Line 1: dimension (1-4)
    - Line 2: size in x-direction (ax)
    - Line 3: size in y-direction (ay) [if dim > 1]
    - Line 4: size in z-direction (az) [if dim > 2]
    - Line 5: size in w-direction (aw) [if dim > 3]
    - Remaining lines: data values (one per line, -1 replaced with threshold)

    Data is stored in Fortran order (column-major).

    Args:
        input_file: Path to input PERSEUS .txt file
        output_file: Path to output .npy file
        verbose: Enable detailed output
    """
    try:
        print(f"Reading PERSEUS format from: {input_file}")

        with open(input_file, 'r') as f:
            lines = f.readlines()

        # Remove empty lines and strip whitespace
        lines = [line.strip() for line in lines if line.strip()]

        if len(lines) < 2:
            raise ValueError("Invalid PERSEUS format: insufficient header lines")

        # Parse header
        line_idx = 0

        # Line 1: dimension
        dim = int(lines[line_idx])
        line_idx += 1

        if dim < 1 or dim > 4:
            raise ValueError(f"Dimension must be 1-4, got {dim}")

        # Line 2: x-dimension size
        ax = int(lines[line_idx])
        line_idx += 1

        # Line 3: y-dimension size (if dim > 1)
        if dim > 1:
            ay = int(lines[line_idx])
            line_idx += 1
        else:
            ay = 1

        # Line 4: z-dimension size (if dim > 2)
        if dim > 2:
            az = int(lines[line_idx])
            line_idx += 1
        else:
            az = 1

        # Line 5: w-dimension size (if dim > 3)
        if dim > 3:
            aw = int(lines[line_idx])
            line_idx += 1
        else:
            aw = 1

        # Calculate expected data size
        expected_size = ax * ay * az * aw
        shape = [ax, ay, az, aw][:dim]  # Only keep relevant dimensions

        if verbose:
            print(f"  Dimension: {dim}")
            print(f"  Shape: {shape}")
            print(f"  Expected data points: {expected_size}")

        # Parse data values
        data_lines = lines[line_idx:]
        if len(data_lines) < expected_size:
            print(f"  Warning: Found {len(data_lines)} data points, expected {expected_size}")

        # Convert data, replacing -1 with a threshold value (using 0 as default)
        threshold = 0.0  # You can make this configurable if needed
        data_values = []

        for i, line in enumerate(data_lines[:expected_size]):
            try:
                value = float(line)
                if value == -1:
                    data_values.append(threshold)
                else:
                    data_values.append(value)
            except ValueError:
                raise ValueError(f"Invalid data value at line {line_idx + i + 1}: '{line}'")

        # Pad with threshold if not enough data
        while len(data_values) < expected_size:
            data_values.append(threshold)

        # Convert to numpy array in Fortran order (column-major)
        data_array = np.array(data_values, dtype=np.float64)

        # Reshape according to PERSEUS format (Fortran order)
        if dim == 1:
            result_array = data_array.reshape((ax,), order='F')
        elif dim == 2:
            result_array = data_array.reshape((ax, ay), order='F')
        elif dim == 3:
            result_array = data_array.reshape((ax, ay, az), order='F')
        else:  # dim == 4
            result_array = data_array.reshape((ax, ay, az, aw), order='F')

        # Display statistics
        print(f"Converted PERSEUS to array: shape={result_array.shape}")
        print(f"  Data range: [{result_array.min():.6f}, {result_array.max():.6f}]")
        print(f"  Mean value: {result_array.mean():.6f}")

        if verbose:
            print(f"  Data type: {result_array.dtype}")
            print(f"  Memory layout: {'Fortran' if result_array.flags['F_CONTIGUOUS'] else 'C'} contiguous")
            unique_vals = np.unique(result_array)
            print(f"  Unique values: {len(unique_vals)} (showing first 10: {unique_vals[:10]})")

        # Save the array
        np.save(output_file, result_array)
        print(f"Saved to: {output_file}")

    except Exception as e:
        print(f"Error converting PERSEUS file: {e}")
        raise

def convert_dipha_persistence_to_numpy(input_file: str, output_file: str, verbose: bool = False) -> None:
    """
    Convert DIPHA persistence diagram to NumPy array.

    DIPHA output/diagram files contain persistence diagrams showing
    the birth and death times of topological features.

    Args:
        input_file: Path to input .output or .diagram file
        output_file: Path to output .npy file
        verbose: Enable detailed output
    """
    try:
        with open(input_file, 'rb') as f:
            data = f.read()

        # Read header: magic number, type, number of persistence pairs
        magic, data_type, num_pairs = struct.unpack_from("qqq", data, 0)

        if verbose:
            print(f"DIPHA Persistence Magic: {magic}, Type: {data_type}, Pairs: {num_pairs}")

        # Read persistence data: (dimension, birth, death) triplets
        offset = 8 * 3  # 8 bytes per long, 3 longs in header
        persistence_data = struct.unpack_from("qdd" * num_pairs, data, offset)

        # Reshape into (n, 3) array: [dimension, birth_time, death_time]
        result_array = np.array(persistence_data).reshape((-1, 3))

        # Count features by dimension
        dims = result_array[:, 0].astype(int)
        dim_counts = {}
        for d in range(int(dims.max()) + 1):
            count = len(result_array[dims == d])
            if count > 0:
                dim_counts[d] = count

        print("Converted persistence diagram:")
        for dim, count in sorted(dim_counts.items()):
            print(f"  {dim}-dimensional features: {count}")
        print(f"Total features: {len(result_array)}")

        if verbose and len(result_array) > 0:
            print(f"First 5 features:")
            for i in range(min(5, len(result_array))):
                dim, birth, death = result_array[i]
                print(f" Dim {int(dim)}: birth={birth:.6f}, death={death:.6f}")

        np.save(output_file, result_array)
        print(f"Saved to: {output_file}")

    except Exception as e:
        print(f"Error converting persistence diagram: {e}")
        sys.exit(1)

def convert_diamorse_to_numpy(input_file: str, output_file: str, verbose: bool = False) -> None:
    """
    Convert DiaMorse output text format to NumPy array.
    Args:
        input_file: Path to input .txt file
        output_file: Path to output .npy file
        verbose: Enable detailed output
    """
    try:
        # Read non-comment lines from file
        text_lines = []
        with open(input_file, 'r') as f:
            for line in f:
                if not line.strip().startswith("#"):
                    text_lines.append(line)

        if not text_lines:
            raise ValueError("No data lines found in DiaMorse file")

        # Parse the data and reorder columns: [dim, birth, death, ...additional_columns]
        raw_data = np.genfromtxt(text_lines)

        # DiaMorse format reordering: take columns [2,0,1,3,4,5,6,7,8]
        if raw_data.shape[1] < 9:
            raise ValueError(f"Expected at least 9 columns in DiaMorse format, got {raw_data.shape[1]}")

        result_array = raw_data[:, [2, 0, 1, 3, 4, 5, 6, 7, 8]]

        # Count features by dimension
        dims = result_array[:, 0].astype(int)
        dim_counts = {}
        for d in range(int(dims.max()) + 1):
            count = len(result_array[dims == d])
            if count > 0:
                dim_counts[d] = count

        print("Converted DiaMorse data:")
        for dim, count in sorted(dim_counts.items()):
            print(f"  {dim}-dimensional features: {count}")
        print(f"  Total features: {len(result_array)}")

        if verbose and len(result_array) > 0:
            print(f"  Data shape: {result_array.shape}")
            print(f"  First 5 features:")
            for i in range(min(5, len(result_array))):
                print(f"    {result_array[i]}")

        np.save(output_file, result_array)
        print(f"Saved to: {output_file}")

    except Exception as e:
        print(f"Error converting DiaMorse file: {e}")
        sys.exit(1)

def detect_file_type(filename: str) -> str:
    """Detect the file type based on filename extension."""
    filename_lower = filename.lower()

    if ".complex" in filename_lower:
        return "dipha_complex"
    elif ".npy" in filename_lower:
        return "numpy"
    elif ".output" in filename_lower or ".diagram" in filename_lower:
        return "dipha_persistence"
    elif ".txt" in filename_lower:
        # Try to distinguish between PERSEUS and DiaMorse formats
        try:
            with open(filename, 'r') as f:
                first_line = f.readline().strip()
                # PERSEUS starts with dimension number (1-4)
                if first_line.isdigit() and 1 <= int(first_line) <= 4:
                    second_line = f.readline().strip()
                    # Second line should also be a positive integer (size)
                    if second_line.isdigit() and int(second_line) > 0:
                        return "perseus"
                # DiaMorse typically starts with comments or has multiple columns
                return "diamorse"
        except:
            return "diamorse"  # Default to diamorse for .txt files
    else:
        return "unknown"

def main() -> None:
    """Main function to handle command line arguments and file conversion."""
    parser = setup_argument_parser()
    args = parser.parse_args()

    # Validate input file exists
    input_path = Path(args.from_file)
    if not input_path.exists():
        print(f"Error: Input file '{args.from_file}' does not exist.")
        sys.exit(1)

    # Detect file types
    input_type = detect_file_type(args.from_file)

    # Create output directory if needed
    output_path = Path(args.to_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Converting: {args.from_file} → {args.to_file}")
    print(f"  Input type: {input_type}")
    print("-" * 50)

    # Perform conversion based on input file type
    try:
        if input_type == "dipha_complex":
            convert_dipha_complex_to_numpy(args.from_file, args.to_file, args.verbose)

        elif input_type == "numpy":
            convert_numpy_to_dipha_complex(args.from_file, args.to_file, args.verbose)

        elif input_type == "dipha_persistence":
            convert_dipha_persistence_to_numpy(args.from_file, args.to_file, args.verbose)

        elif input_type == "perseus":
            convert_perseus_to_numpy(args.from_file, args.to_file, args.verbose)

        elif input_type == "diamorse":
            convert_diamorse_to_numpy(args.from_file, args.to_file, args.verbose)

        else:
            print(f"Error: Unsupported file type for '{args.from_file}'")
            print("   Supported formats: .complex, .npy, .output, .diagram, .txt")
            sys.exit(1)

        print("-" * 50)
        print("Conversion completed successfully!")

    except KeyboardInterrupt:
        print("\n  Conversion interrupted by user.")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
