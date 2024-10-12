from abc import ABC





class RetrievableSettings(ABC):

    def _class_name(self) -> str:
        class_name = self.__class__.__qualname__
        module_name = self.__class__.__module__
        result = f'{module_name}.{class_name}'
        if module_name == 'builtins': # Ignore 'builtins.' module name
            result = class_name
        i = result.find('.') # Remove 'evodesign.' at the beginning
        j = result.rfind('.') # Remove repeated class name at the end
        return result[i + 1 : j]
  


    def settings(self) -> dict:
        # returns all the public instance and class attributes as a dictionary
        settings_dict = {}
        for key, value in self.__dict__.items():
            if isinstance(value, RetrievableSettings):
                # apply the function recursively
                settings_dict[key] = value.settings()
            elif not key.startswith('_'):
                settings_dict[key] = value
        return { self._class_name(): settings_dict }
