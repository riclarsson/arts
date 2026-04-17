import pyarts3 as pyarts

ws = pyarts.Workspace()

numeric1 = pyarts.arts.Numeric(2.0)
numeric2 = pyarts.arts.Numeric(4.0)

# Test that WriteXML and ReadXML work for Numeric (value-type)
ws.WriteXML(filename="test.xml", input=numeric1)
ws.ReadXML(numeric2, filename="test.xml")
assert numeric1 == numeric2
assert numeric2 == 2.0

aov1 = pyarts.arts.ArrayOfAscendingGrid([[1, 2, 3], [1, 2]])
aov2 = pyarts.arts.ArrayOfAscendingGrid([[4, 5, 6], [4, 5]])

# Test that WriteXML and ReadXML work for ArrayOfAscendingGrid (reference-type)
ws.WriteXML(filename="test.xml", input=aov1)
ws.ReadXML(aov2, filename="test.xml")

assert aov1 == aov2
assert aov2 == [[1, 2, 3], [1, 2]]
