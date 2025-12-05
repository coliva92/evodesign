import importlib
import copy
import json


def load(settings_path: str):
    with open(settings_path, "rt", encoding="utf-8") as json_file:
        settings = json.load(json_file)
    return parse(settings)


def parse(settings: dict):
    class_name = list(settings.keys())[0]
    imported_module = importlib.import_module(f"evodesign.{class_name}")
    actual_class = getattr(imported_module, class_name.split(".")[-1])
    params = copy.deepcopy(settings[class_name])
    for key, item in params.items():
        if type(item) == dict:
            params[key] = parse(item)
        if type(item) == list and type(item[0]) == dict:
            params[key] = [parse(s) for s in item] if type(item[0]) == dict else item
    return actual_class(**params)
