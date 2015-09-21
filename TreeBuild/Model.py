'''
    Keep all consistent variable
'''


#SMILE_COLUMNNAME = "Smiles"
SMILE_COLUMNNAME = "Canonical_Smiles"
#SMILE_COLUMNNAME = "Smiles"

#FILE_FORMAT = './Data/%Y-%m-%d-%Hh-%Mm'
FILE_FORMAT = './Data/%Y-%m-%d-%Hh-%Mm-%Ss'

# image directory
IMG_DIR = "./Data/Image/"

# Potency
POTENCY = "pIC50"

# interesting column
#INTEREST_COLUMN = ["IC50", "orig_id"]
INTEREST_COLUMN = ["pIC50", "size", "group", "orig_id"]
ACTIVITY_COLUMN = ["pIC50"]

# identifier
IDENTIFIER = "ligandid"
