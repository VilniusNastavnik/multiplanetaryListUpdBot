import unittest

import multiplanetaryListUpdBot as multiplanetaryListUpdBot

class mainCheck(unittest.TestCase):

    def setUp(self):
        pass

    def test_deg_to_hms(self):

        ra = multiplanetaryListUpdBot.deg_to_hms("164.120833343","RA")
        dec = multiplanetaryListUpdBot.deg_to_hms("7.014444456","DEC")

        self.assertTrue(ra == "10|56|29",ra)
        self.assertTrue(dec == "7|0|52",dec)

    def test_hms_to_numb(self):
        
        n = multiplanetaryListUpdBot.hms_to_numb("7|1|52")
        self.assertTrue(n == 25312,n)

    def test_getCoordFromSimbad(self):

        dist, ra, dec = multiplanetaryListUpdBot.getCoordFromSimbad("24 Sex")
        self.assertTrue(ra == "10|23|28",ra)
        self.assertTrue(dec == "0|54|8",dec)
        self.assertTrue(dist == 239.1,dist)

if __name__ == "__main__":
    unittest.main()