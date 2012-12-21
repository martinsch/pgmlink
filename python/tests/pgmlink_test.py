import sys
sys.path.append("../.")

import unittest as ut
import pgmlink

def mk_traxel(x,y,z,id):
    t = pgmlink.Traxel()
    t.ID = id
    t.add_feature_array("com", 3)
    t.set_feature_value("com", 0, x)
    t.set_feature_value("com", 1, y)
    t.set_feature_value("com", 2, z)
    return t

class Test_Traxels( ut.TestCase ):
    def runTest( self ):
        t = mk_traxel(34,23,12,77)
        self.assertEqual(t.ID, 77)
        self.assertEqual(len(t.features), 1)
        self.assertEqual(t.get_feature_value("com", 0), 34)
        self.assertEqual(t.get_feature_value("com", 1), 23)
        self.assertEqual(t.get_feature_value("com", 2), 12)

if __name__=="__main__":
    ut.main()
