import os
import re
from urllib.parse import urljoin

def convert_to_datasets_path(local_path: str, base_url: str = "https://genome.med.nyu.edu/") -> str:
    """
    Converts a local file system path to a URL suitable for datasets.json.

    Assumes a structure like /gpfs/data/{lab_name}/external/{rest_of_path}
    and converts it to {base_url}/external/{lab_name}/{rest_of_path}/

    Args:
        local_path: The absolute local file system path.
        base_url: The base URL for the web server.

    Returns:
        The converted URL path.
    """
    # Normalize the path to remove trailing slashes
    local_path = os.path.normpath(local_path)

    # Regex to capture the lab name and the rest of the path after 'external'
    # It looks for '/{lab_name}/external/' pattern
    match = re.search(r'/([^/]+)/external/(.*)', local_path)

    if not match:
        raise ValueError(f"Path '{local_path}' does not match the expected pattern '.../{{lab_name}}/external/{{rest_of_path}}'")

    lab_name = match.group(1)
    rest_of_path = match.group(2)

    # Construct the relative path part for the URL
    # Ensure it ends with a slash
    relative_url_path = os.path.join('external', lab_name, rest_of_path).replace(os.sep, '/') + '/'

    # Join the base URL with the constructed relative path
    final_url = urljoin(base_url, relative_url_path)

    return final_url

# Example usage:
if __name__ == '__main__':
    test_path = "/gpfs/data/imielinskilab/external/mskiweb/diders01/gos_hmf_test"
    converted_url = convert_to_datasets_path(test_path)
    print(f"Original path: {test_path}")
    print(f"Converted URL: {converted_url}")
    # Expected output: https://genome.med.nyu.edu/external/imielinskilab/mskiweb/diders01/gos_hmf_test/

    test_path_trailing_slash = "/gpfs/data/imielinskilab/external/mskiweb/diders01/gos_hmf_test/"
    converted_url_trailing = convert_to_datasets_path(test_path_trailing_slash)
    print(f"Original path: {test_path_trailing_slash}")
    print(f"Converted URL: {converted_url_trailing}")
    # Expected output: https://genome.med.nyu.edu/external/imielinskilab/mskiweb/diders01/gos_hmf_test/

    # Test with a different lab
    test_path_other_lab = "/gpfs/data/otherlab/external/projectx/data"
    converted_url_other_lab = convert_to_datasets_path(test_path_other_lab)
    print(f"Original path: {test_path_other_lab}")
    print(f"Converted URL: {converted_url_other_lab}")
    # Expected output: https://genome.med.nyu.edu/external/otherlab/projectx/data/

    # Test invalid path
    try:
        invalid_path = "/some/other/path/structure"
        convert_to_datasets_path(invalid_path)
    except ValueError as e:
        print(f"Caught expected error for invalid path: {e}")
