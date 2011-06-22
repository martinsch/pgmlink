import sys
sys.path.append("../.")

import unittest as ut
import ctracking as ct

def mk_tracklet(x,y,z,id):
    t = ct.Tracklet()
    t.ID = id
    t.add_feature_array("com", 3)
    t.set_feature_value("com", 0, x)
    t.set_feature_value("com", 1, y)
    t.set_feature_value("com", 2, z)
    return t

class Test_Tracklet( ut.TestCase ):
    def runTest( self ):
        t = mk_tracklet(34,23,12,77)
        self.assertEqual(t.ID, 77)
        self.assertEqual(len(t.features), 1)
        self.assertEqual(t.get_feature_value("com", 0), 34)
        self.assertEqual(t.get_feature_value("com", 1), 23)
        self.assertEqual(t.get_feature_value("com", 2), 12)

class Test_AtomicEvents( ut.TestCase ):
    def setUp( self ):
        self.prev = ct.Tracklets()
        self.prev.add_tracklet(mk_tracklet(0,0,0,4))
        self.curr = ct.Tracklets()
        self.curr.add_tracklet(mk_tracklet(1,0,0,12))
    
    def test_disappearance( self ):
        m = ct.FixedCostTracking(10,10,10,10,50)
        events = m(self.prev,ct.Tracklets())
        self.assertEqual(len(events), 1)
        self.assertEqual(len(events[0].tracklet_ids), 1)
        self.assertEqual(events[0].tracklet_ids[0], 4)
        self.assertEqual(events[0].type, ct.EventType.Disappearance)

    def test_move( self ):
        m = ct.FixedCostTracking(10,10,1000,1000,50)
        events = m(self.prev, self.curr)
        self.assertEqual(len(events), 1)
        self.assertEqual(len(events[0].tracklet_ids), 2)
        self.assertEqual(events[0].tracklet_ids[0], 4)
        self.assertEqual(events[0].tracklet_ids[1], 12)
        self.assertEqual(events[0].energy, 10)

if __name__=="__main__":
    ut.main()
