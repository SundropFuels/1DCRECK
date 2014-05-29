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

    def __add__(self, other):
        """Adds two unitvals"""
        return UnitVal(self.value+other.get_value(self.units), self.units)

    def __radd__(self, other):
        return self.__add__(self, other)

    def __sub__(self, other):
        """Subtracts two unitvals"""
        return UnitVal(self.value - other.get_value(self.units), self.units)

    def __rsub__(self, other):
        return self.__sub__(self, other)


    def __mul__(self, other):
        """Multiplies a unit value with another unit value or a scalar"""
        if isinstance(other, UnitVal):
            return UnitVal(self.value*other.get_value(self.units), self.units)
        else:
            return UnitVal(self.value*other, self.units)

    def __rmul__(self, other):
        return self.__mul__(self, other)


    def __div__(self, other):
        """Divides a unit value with another unit value or a scalar"""

       if isinstance(other, UnitVal):
           if other.value == 0:
               raise ValueError, "Divide by zero"
           else:
               return UnitVal(self.value*other.get_value(self.units), self.units)
       else:
           if other.value == 0:
               raise ValueError, "Divide by zero"
           else:
               return UnitVal(self.value*other, self.units)

    def __rdiv__(self, other):
        return self.__div__(self, other)
