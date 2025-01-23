import ast

# Function to extract variable and function names from Python code
def extract_names_from_file(file_path):
    with open(file_path, 'r') as file:
        code = file.read()

    # Parse the code into an abstract syntax tree (AST)
    tree = ast.parse(code)

    # Define sets to hold variable and function names
    var_names = set()
    func_names = set()

    # Walk through all nodes in the AST
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):  # Check for assignment statements
            for target in node.targets:
                if isinstance(target, ast.Name):  # Check if the target is a variable
                    var_names.add(target.id)
        elif isinstance(node, ast.FunctionDef):  # Check for function definitions
            func_names.add(node.name)

    # Print the variable names and function names separated by commas
    print("Variables:", ", ".join(sorted(var_names)))
    print("Functions:", ", ".join(sorted(func_names)))

# Call the function with the path to your Python file
file_path = 'Beam_intensity.py'  # Replace with the actual file path
extract_names_from_file(file_path)
