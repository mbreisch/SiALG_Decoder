import sys 
import uproot as ur4


"""
Read a root file using uproot4 and extract data from specified trees and branches.

Args:
    path (str): The path to the root file.
    treenames (list): A list of tree names to read.
    branchnames (list): A list of branch names to extract data from.

Returns:
    dict: A dictionary containing the structured data.

"""
def read_root_file(path, treenames, branchnames):
    try:
        file = ur4.open(path)
        data = {}
        for treename in treenames:
            tree = file[treename]
            data[treename] = {}
            for branchname in branchnames:
                branch = tree[branchname]
                data[treename][branchname] = branch.array()[0]
        return data
    except FileNotFoundError:
        print(f"File '{path}' not found.")
        return False
    except KeyError:
        #print(f"Tree or branch not found in file '{path}'.")
        return False
    except Exception as e:
        print(f"An error occurred while reading the root file: {str(e)}")
        return False


"""
Read a root file dictionary and prints trees and branches.

Args:
    roottree (dictionary): The compiled rootfile dictionary from read_root_file

Returns:
    None

"""
def PrintRootTrees(roottree):
    if roottree is not None:
        for treename, branches in roottree.items():
            print(f"Tree: {treename}")
            for branchname, data in branches.items():
                print(f"\tBranch: {branchname}")
                print("\t     ",data)


"""
Reads a list of path strings from a file.

Parameters:
    file_path (str): The path to the file containing the list of path strings.

Returns:
    list: A list of path strings read from the file.

Raises:
    FileNotFoundError: If the specified file cannot be found.
    IOError: If there is an error while reading the file.

"""
def load_file_list():
    file_path = "./path_list_jupyter"
    try:
        with open(file_path, 'r') as file:
            path_list = file.read().splitlines()
    except FileNotFoundError:
        raise FileNotFoundError(f"File '{file_path}' not found.")
    except IOError as e:
        raise IOError(f"Error while reading the file '{file_path}': {str(e)}")

    return path_list


"""
Read a root file using uproot4 and extract data from specified trees and branches.

Args:
    path (str): The path to the root file.
    treenames (list): A list of tree names to read.
    branchnames (list): A list of branch names to extract data from.

Returns:
    dict: A dictionary containing the structured data.

"""
def load_events(path, treenames, branchnames):
    try:
        file = ur4.open(path)
        data = {}
        for treename in treenames:
            tree = file[treename]
            data[treename] = {}
            for branchname in branchnames:
                branch = tree[branchname]
                data[treename][branchname] = branch.array()[0]
        return data
    except FileNotFoundError:
        print(f"File '{path}' not found.")
        return False
    except KeyError:
        #print(f"Tree or branch not found in file '{path}'.")
        return False
    except Exception as e:
        print(f"An error occurred while reading the root file: {str(e)}")
        return False

