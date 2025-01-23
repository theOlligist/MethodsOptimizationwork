import sys
import re


def make_seq_dictionary(file_path):
    """
    Calculate the total length of valid DNA sequences in the given file.
    
    :param file_path: Path to the file containing DNA sequences.
    :return: Total length of valid sequences as an integer.
    """
    lengths = []
    try:
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()  # Remove leading/trailing whitespace
                # Check if the line contains only valid nucleotide characters
                if re.fullmatch(r"[ATCG]+", line):
                    lengths.append(len(line))
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

    return sum(lengths)


def extract_total_bp(log_file):
    """
    Extract the total number of base pairs (bp) from a MEGAHIT assembly log.

    :param log_file: Path to the MEGAHIT assembly log file.
    :return: Total number of base pairs as an integer, or None if not found.
    """
    try:
        with open(log_file, "r") as file:
            for line in file:
                match = re.search(r"total\s+([\d,]+)\s+bp,\s+min", line)
                if match:
                    return int(match.group(1).replace(",", ""))
    except FileNotFoundError:
        print(f"Error: The file '{log_file}' was not found.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

    return None


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python combined_script.py <log_file_path> <sequence_file_path>")
        sys.exit(1)

    log_file = sys.argv[1]
    sequence_file = sys.argv[2]

    # Extract total_bp from the log file
    total_bp = extract_total_bp(log_file)
    if total_bp is None:
        print("Unable to extract total base pairs from the log file.")
        sys.exit(1)

    # Calculate total_length from the sequence file
    total_length = make_seq_dictionary(sequence_file)
    if total_length is None:
        print("Unable to calculate total length of valid sequences.")
        sys.exit(1)

    # Calculate and display the ratio
    if total_length == 0:
        print("Error: Total length of sequences is zero, cannot compute ratio.")
    else:
        ratio = total_bp / total_length
        print(f"Total BP Assembly: {total_bp}")
        print(f"Total BP Reads: {total_length}")
        print(f"Ratio (total_bp / total_length): {ratio:.4f}")