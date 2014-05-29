import unitConversion as uc



class UnitVal:
    def __init__(self, value = None, units = None):
        self.value = value
        self.units = units

    def set_unit_class(self):
        """Will find out to what class of physical measurements (e.g. length, time, pressure, temperature, etc.) the units belong"""
        pass

    def __eq__(self, other):
        return isinstance(other, self.__class__) and self.value == other.value and self.units == other.units

    def __ne__(self, other):
        return not self.__eq__(other)

    def get_value(self, unit):
        """Returns the value in the requested unit"""
        conv = uc.UnitConverter()
        return conv.convert_units(self.value, self.units, unit)   #This needs error handling
