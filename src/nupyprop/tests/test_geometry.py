from nupyprop import geometry
from pytest import approx

def test_find_interface():
   assert geometry.find_interface(1) == approx([1.01517171420510, 358.984828285795], rel=1e-4)
   assert geometry.find_interface(0) == approx([0.0, 360.0], rel=1e-4)
