import sys

root_system = sys.argv[1]
representation = int(sys.argv[2])

R = RootSystem(root_system)
A = R.ambient_space()

v_0 = A.fundamental_weight(representation)
weyl_orbit = v_0.orbit()

data = {
    "root_system": root_system,
    "representation": representation,
    "weights": weyl_orbit,
    "roots": A.roots(),
}

print data
