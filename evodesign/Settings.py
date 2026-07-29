import copy
import importlib
import json
import os

from .AlgorithmFactory import AlgorithmFactory


def read_json(settings_path: str) -> dict:
    with open(os.path.abspath(settings_path), "rt", encoding="utf-8") as json_file:
        settings = json.load(json_file)
    return settings


def parse(settings: dict) -> AlgorithmFactory:
    # obtain the full namespace resolution string and import the class
    class_namespace = list(settings.keys())[0]
    imported_module = importlib.import_module(f"evodesign.{class_namespace}")

    # obtain the actual class definition and its constructor arguments
    actual_class = getattr(imported_module, class_namespace.split(".")[-1])
    constructor_args = copy.deepcopy(settings[class_namespace])

    # parse each constructor argument if it's a nested class or a class list
    for key, item in constructor_args.items():
        if type(item) == dict:
            constructor_args[key] = parse(item)
        if type(item) == list and type(item[0]) == dict:
            constructor_args[key] = (
                [parse(s) for s in item] if type(item[0]) == dict else item
            )

    # return the instance of the root class
    return actual_class(**constructor_args)


def load(settings_path: str) -> AlgorithmFactory:
    settings = read_json(settings_path)
    return parse(settings)


def get_namespace(obj):
    """Extracts the full namespace resolution string of an object."""
    class_name = obj.__class__.__qualname__
    module_name = obj.__class__.__module__
    result = f"{module_name}.{class_name}"

    if module_name != "builtins":  # ignore 'builtins.' module name
        i = result.find(".")  # remove 'evodesign.' at the beginning
        j = result.rfind(".")  # remove repeated class name at the end
        if i != -1 and j != -1 and i < j:
            result = result[i + 1 : j]
    return result


class _SettingsEncoder(json.JSONEncoder):

    def default(self, obj):
        # dump the public attributes and constructor arguments
        # basic primitive types are handled automatically by (JSONEncoder)
        settings_dict = {
            key: value
            for key, value in obj.__dict__.items()
            if not key.startswith("_")  # ignore private attributes
        }
        return {get_namespace(obj): settings_dict}


def dump(obj) -> dict:
    return json.loads(json.dumps(obj, cls=_SettingsEncoder))
