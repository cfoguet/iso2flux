import math

def round_up(number,positions):
    exponent=pow(10,positions)
    new_number=math.ceil(number*exponent)/exponent
    """if new_number==number:
       new_number=number+1.0/exponent"""
    return new_number


def round_down(number,positions):
    if number==0.0:
       return 0
    exponent=pow(10,positions)
    return math.floor(number*exponent-0.0001)/exponent
    """if new_number==number:
       new_number=number-1.0/exponent"""
    return new_number
