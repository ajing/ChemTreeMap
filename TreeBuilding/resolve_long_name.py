'''
 resolve long name problem
'''

from Dot2JSON import Dot2JSON, Root2JSON

root = Dot2JSON("test_sfdp.gv")
Root2JSON(root, "test.json")
