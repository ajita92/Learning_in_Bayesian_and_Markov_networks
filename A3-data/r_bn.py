from bayesian.factor_graph import *
from bayesian.bbn import *

dictionary_GOAL_2 = {'true': 1.0, 'false': 1.0}

def f_GOAL_2(GOAL_2):
    return dictionary_GOAL_2[GOAL_2]

dictionary_SNode_3 = {'true': 1.0, 'false': 1.0}

def f_SNode_3(SNode_3):
    return dictionary_SNode_3[SNode_3]

dictionary_SNode_4 = {'true': 1.0, 'false': 1.0}

def f_SNode_4(SNode_4):
    return dictionary_SNode_4[SNode_4]

dictionary_SNode_5 = {'true': 1.0, 'false': 1.0}

def f_SNode_5(SNode_5):
    return dictionary_SNode_5[SNode_5]

dictionary_SNode_6 = {'true': 1.0, 'false': 1.0}

def f_SNode_6(SNode_6):
    return dictionary_SNode_6[SNode_6]

dictionary_SNode_7 = {'true': 1.0, 'false': 1.0}

def f_SNode_7(SNode_7):
    return dictionary_SNode_7[SNode_7]

dictionary_DISPLACEM0 = {'true': 1.0, 'false': 1.0}

def f_DISPLACEM0(DISPLACEM0):
    return dictionary_DISPLACEM0[DISPLACEM0]

dictionary_RApp1 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_RApp1(DISPLACEM0, SNode_3, RApp1):
    return dictionary_RApp1[(DISPLACEM0, SNode_3, RApp1)]

dictionary_GIVEN_1 = {'true': 1.0, 'false': 1.0}

def f_GIVEN_1(GIVEN_1):
    return dictionary_GIVEN_1[GIVEN_1]

dictionary_RApp2 = {('false', 'true'): 1.0, ('true', 'false'): 1.0, ('false', 'false'): 1.0, ('true', 'true'): 1.0}
def f_RApp2(GIVEN_1, RApp2):
    return dictionary_RApp2[(GIVEN_1, RApp2)]

dictionary_SNode_8 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_8(RApp1, RApp2, SNode_8):
    return dictionary_SNode_8[(RApp1, RApp2, SNode_8)]

dictionary_SNode_9 = {'true': 1.0, 'false': 1.0}

def f_SNode_9(SNode_9):
    return dictionary_SNode_9[SNode_9]

dictionary_SNode_10 = {'true': 1.0, 'false': 1.0}

def f_SNode_10(SNode_10):
    return dictionary_SNode_10[SNode_10]

dictionary_SNode_11 = {'true': 1.0, 'false': 1.0}

def f_SNode_11(SNode_11):
    return dictionary_SNode_11[SNode_11]

dictionary_SNode_12 = {'true': 1.0, 'false': 1.0}

def f_SNode_12(SNode_12):
    return dictionary_SNode_12[SNode_12]

dictionary_SNode_13 = {'true': 1.0, 'false': 1.0}

def f_SNode_13(SNode_13):
    return dictionary_SNode_13[SNode_13]

dictionary_SNode_14 = {'true': 1.0, 'false': 1.0}

def f_SNode_14(SNode_14):
    return dictionary_SNode_14[SNode_14]

dictionary_SNode_15 = {'true': 1.0, 'false': 1.0}

def f_SNode_15(SNode_15):
    return dictionary_SNode_15[SNode_15]

dictionary_SNode_16 = {'true': 1.0, 'false': 1.0}

def f_SNode_16(SNode_16):
    return dictionary_SNode_16[SNode_16]

dictionary_SNode_17 = {'true': 1.0, 'false': 1.0}

def f_SNode_17(SNode_17):
    return dictionary_SNode_17[SNode_17]

dictionary_SNode_18 = {'true': 1.0, 'false': 1.0}

def f_SNode_18(SNode_18):
    return dictionary_SNode_18[SNode_18]

dictionary_SNode_19 = {'true': 1.0, 'false': 1.0}

def f_SNode_19(SNode_19):
    return dictionary_SNode_19[SNode_19]

dictionary_NEED1 = {'true': 1.0, 'false': 1.0}

def f_NEED1(NEED1):
    return dictionary_NEED1[NEED1]

dictionary_SNode_20 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_20(SNode_16, NEED1, SNode_20):
    return dictionary_SNode_20[(SNode_16, NEED1, SNode_20)]

dictionary_GRAV2 = {'true': 1.0, 'false': 1.0}

def f_GRAV2(GRAV2):
    return dictionary_GRAV2[GRAV2]

dictionary_SNode_21 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_21(SNode_20, GRAV2, SNode_21):
    return dictionary_SNode_21[(SNode_20, GRAV2, SNode_21)]

dictionary_VALUE3 = {'true': 1.0, 'false': 1.0}

def f_VALUE3(VALUE3):
    return dictionary_VALUE3[VALUE3]

dictionary_SNode_24 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_24(SNode_21, VALUE3, SNode_24):
    return dictionary_SNode_24[(SNode_21, VALUE3, SNode_24)]

dictionary_SLIDING4 = {'true': 1.0, 'false': 1.0}

def f_SLIDING4(SLIDING4):
    return dictionary_SLIDING4[SLIDING4]

dictionary_SNode_25 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_25(SNode_15, SLIDING4, SNode_25):
    return dictionary_SNode_25[(SNode_15, SLIDING4, SNode_25)]

dictionary_CONSTANT5 = {'true': 1.0, 'false': 1.0}

def f_CONSTANT5(CONSTANT5):
    return dictionary_CONSTANT5[CONSTANT5]

dictionary_SNode_26 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_26(SNode_11, CONSTANT5, SNode_26):
    return dictionary_SNode_26[(SNode_11, CONSTANT5, SNode_26)]

dictionary_KNOWN6 = {'true': 1.0, 'false': 1.0}

def f_KNOWN6(KNOWN6):
    return dictionary_KNOWN6[KNOWN6]

dictionary_VELOCITY7 = {'true': 1.0, 'false': 1.0}

def f_VELOCITY7(VELOCITY7):
    return dictionary_VELOCITY7[VELOCITY7]

dictionary_SNode_47 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_47(SNode_3, VELOCITY7, SNode_47):
    return dictionary_SNode_47[(SNode_3, VELOCITY7, SNode_47)]

dictionary_RApp3 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_RApp3(KNOWN6, SNode_26, SNode_47, RApp3):
    return dictionary_RApp3[(KNOWN6, SNode_26, SNode_47, RApp3)]

dictionary_KNOWN8 = {'true': 1.0, 'false': 1.0}

def f_KNOWN8(KNOWN8):
    return dictionary_KNOWN8[KNOWN8]

dictionary_RApp4 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_RApp4(KNOWN8, SNode_11, RApp4):
    return dictionary_RApp4[(KNOWN8, SNode_11, RApp4)]

dictionary_SNode_27 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_27(RApp3, RApp4, SNode_27):
    return dictionary_SNode_27[(RApp3, RApp4, SNode_27)]

dictionary_COMPO16 = {'true': 1.0, 'false': 1.0}

def f_COMPO16(COMPO16):
    return dictionary_COMPO16[COMPO16]

dictionary_GOAL_48 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_GOAL_48(GOAL_2, COMPO16, GOAL_48):
    return dictionary_GOAL_48[(GOAL_2, COMPO16, GOAL_48)]

dictionary_TRY12 = {'true': 1.0, 'false': 1.0}

def f_TRY12(TRY12):
    return dictionary_TRY12[TRY12]

dictionary_TRY11 = {('false', 'true'): 1.0, ('true', 'false'): 1.0, ('false', 'false'): 1.0, ('true', 'true'): 1.0}
def f_TRY11(TRY12, TRY11):
    return dictionary_TRY11[(TRY12, TRY11)]

dictionary_GOAL_49 = {('true', 'true', 'false', 'true', 'true'): 1.0, ('true', 'true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true', 'true'): 1.0, ('true', 'true', 'false', 'false', 'true'): 1.0, ('false', 'true', 'true', 'true', 'true'): 1.0, ('false', 'true', 'true', 'true', 'false'): 1.0, ('false', 'true', 'false', 'false', 'true'): 1.0, ('false', 'false', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'false'): 1.0, ('false', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'false', 'true', 'false'): 1.0, ('true', 'true', 'true', 'true', 'true'): 1.0, ('false', 'true', 'false', 'true', 'false'): 1.0, ('true', 'false', 'true', 'false', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true'): 1.0, ('true', 'false', 'false', 'false', 'false'): 1.0, ('true', 'true', 'false', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'false', 'false', 'true'): 1.0, ('true', 'false', 'false', 'false', 'true'): 1.0, ('false', 'false', 'true', 'true', 'true'): 1.0, ('true', 'true', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true', 'true'): 1.0, ('true', 'true', 'true', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'false', 'false', 'false'): 1.0}
def f_GOAL_49(SNode_5, SNode_6, GOAL_48, TRY11, GOAL_49):
    return dictionary_GOAL_49[(SNode_5, SNode_6, GOAL_48, TRY11, GOAL_49)]

dictionary_CHOOSE19 = {'true': 1.0, 'false': 1.0}

def f_CHOOSE19(CHOOSE19):
    return dictionary_CHOOSE19[CHOOSE19]

dictionary_GOAL_50 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_GOAL_50(GOAL_49, CHOOSE19, GOAL_50):
    return dictionary_GOAL_50[(GOAL_49, CHOOSE19, GOAL_50)]

dictionary_SYSTEM18 = {'true': 1.0, 'false': 1.0}

def f_SYSTEM18(SYSTEM18):
    return dictionary_SYSTEM18[SYSTEM18]

dictionary_SNode_51 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_SNode_51(SNode_17, GOAL_50, SYSTEM18, SNode_51):
    return dictionary_SNode_51[(SNode_17, GOAL_50, SYSTEM18, SNode_51)]

dictionary_KINEMATI17 = {'true': 1.0, 'false': 1.0}

def f_KINEMATI17(KINEMATI17):
    return dictionary_KINEMATI17[KINEMATI17]

dictionary_SNode_52 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_52(SNode_51, KINEMATI17, SNode_52):
    return dictionary_SNode_52[(SNode_51, KINEMATI17, SNode_52)]

dictionary_IDENTIFY10 = {'true': 1.0, 'false': 1.0}

def f_IDENTIFY10(IDENTIFY10):
    return dictionary_IDENTIFY10[IDENTIFY10]

dictionary_GOAL_53 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_GOAL_53(GOAL_49, SNode_52, IDENTIFY10, GOAL_53):
    return dictionary_GOAL_53[(GOAL_49, SNode_52, IDENTIFY10, GOAL_53)]

dictionary_IDENTIFY9 = {'true': 1.0, 'false': 1.0}

def f_IDENTIFY9(IDENTIFY9):
    return dictionary_IDENTIFY9[IDENTIFY9]

dictionary_SNode_28 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_SNode_28(SNode_27, GOAL_53, IDENTIFY9, SNode_28):
    return dictionary_SNode_28[(SNode_27, GOAL_53, IDENTIFY9, SNode_28)]

dictionary_TRY13 = {('false', 'true'): 1.0, ('true', 'false'): 1.0, ('false', 'false'): 1.0, ('true', 'true'): 1.0}
def f_TRY13(TRY12, TRY13):
    return dictionary_TRY13[(TRY12, TRY13)]

dictionary_TRY14 = {('false', 'true'): 1.0, ('true', 'false'): 1.0, ('false', 'false'): 1.0, ('true', 'true'): 1.0}
def f_TRY14(TRY12, TRY14):
    return dictionary_TRY14[(TRY12, TRY14)]

dictionary_TRY15 = {('false', 'true'): 1.0, ('true', 'false'): 1.0, ('false', 'false'): 1.0, ('true', 'true'): 1.0}
def f_TRY15(TRY12, TRY15):
    return dictionary_TRY15[(TRY12, TRY15)]

dictionary_VAR20 = {'true': 1.0, 'false': 1.0}

def f_VAR20(VAR20):
    return dictionary_VAR20[VAR20]

dictionary_SNode_29 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_29(SNode_28, VAR20, SNode_29):
    return dictionary_SNode_29[(SNode_28, VAR20, SNode_29)]

dictionary_SNode_31 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_31(SNode_29, VALUE3, SNode_31):
    return dictionary_SNode_31[(SNode_29, VALUE3, SNode_31)]

dictionary_GIVEN21 = {'true': 1.0, 'false': 1.0}

def f_GIVEN21(GIVEN21):
    return dictionary_GIVEN21[GIVEN21]

dictionary_SNode_33 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_33(SNode_10, GIVEN21, SNode_33):
    return dictionary_SNode_33[(SNode_10, GIVEN21, SNode_33)]

dictionary_SNode_34 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_34(SNode_10, CONSTANT5, SNode_34):
    return dictionary_SNode_34[(SNode_10, CONSTANT5, SNode_34)]

dictionary_VECTOR27 = {'true': 1.0, 'false': 1.0}

def f_VECTOR27(VECTOR27):
    return dictionary_VECTOR27[VECTOR27]

dictionary_APPLY32 = {'true': 1.0, 'false': 1.0}

def f_APPLY32(APPLY32):
    return dictionary_APPLY32[APPLY32]

dictionary_GOAL_56 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_GOAL_56(GOAL_49, SNode_52, APPLY32, GOAL_56):
    return dictionary_GOAL_56[(GOAL_49, SNode_52, APPLY32, GOAL_56)]

dictionary_CHOOSE35 = {'true': 1.0, 'false': 1.0}

def f_CHOOSE35(CHOOSE35):
    return dictionary_CHOOSE35[CHOOSE35]

dictionary_GOAL_57 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_GOAL_57(GOAL_56, CHOOSE35, GOAL_57):
    return dictionary_GOAL_57[(GOAL_56, CHOOSE35, GOAL_57)]

dictionary_MAXIMIZE34 = {'true': 1.0, 'false': 1.0}

def f_MAXIMIZE34(MAXIMIZE34):
    return dictionary_MAXIMIZE34[MAXIMIZE34]

dictionary_SNode_59 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_SNode_59(SNode_7, GOAL_57, MAXIMIZE34, SNode_59):
    return dictionary_SNode_59[(SNode_7, GOAL_57, MAXIMIZE34, SNode_59)]

dictionary_AXIS33 = {'true': 1.0, 'false': 1.0}

def f_AXIS33(AXIS33):
    return dictionary_AXIS33[AXIS33]

dictionary_SNode_60 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_60(SNode_59, AXIS33, SNode_60):
    return dictionary_SNode_60[(SNode_59, AXIS33, SNode_60)]

dictionary_WRITE31 = {'true': 1.0, 'false': 1.0}

def f_WRITE31(WRITE31):
    return dictionary_WRITE31[WRITE31]

dictionary_GOAL_61 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_GOAL_61(GOAL_56, SNode_60, WRITE31, GOAL_61):
    return dictionary_GOAL_61[(GOAL_56, SNode_60, WRITE31, GOAL_61)]

dictionary_WRITE30 = {'true': 1.0, 'false': 1.0}

def f_WRITE30(WRITE30):
    return dictionary_WRITE30[WRITE30]

dictionary_GOAL_62 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_GOAL_62(GOAL_61, WRITE30, GOAL_62):
    return dictionary_GOAL_62[(GOAL_61, WRITE30, GOAL_62)]

dictionary_RESOLVE37 = {'true': 1.0, 'false': 1.0}

def f_RESOLVE37(RESOLVE37):
    return dictionary_RESOLVE37[RESOLVE37]

dictionary_GOAL_63 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_GOAL_63(SNode_28, GOAL_62, RESOLVE37, GOAL_63):
    return dictionary_GOAL_63[(SNode_28, GOAL_62, RESOLVE37, GOAL_63)]

dictionary_NEED36 = {'true': 1.0, 'false': 1.0}

def f_NEED36(NEED36):
    return dictionary_NEED36[NEED36]

dictionary_SNode_64 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_64(GOAL_63, NEED36, SNode_64):
    return dictionary_SNode_64[(GOAL_63, NEED36, SNode_64)]

dictionary_SNode_41 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_41(SNode_9, CONSTANT5, SNode_41):
    return dictionary_SNode_41[(SNode_9, CONSTANT5, SNode_41)]

dictionary_SNode_42 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_SNode_42(SNode_8, SNode_41, KNOWN6, SNode_42):
    return dictionary_SNode_42[(SNode_8, SNode_41, KNOWN6, SNode_42)]

dictionary_IDENTIFY39 = {'true': 1.0, 'false': 1.0}

def f_IDENTIFY39(IDENTIFY39):
    return dictionary_IDENTIFY39[IDENTIFY39]

dictionary_SNode_43 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_SNode_43(SNode_42, GOAL_53, IDENTIFY39, SNode_43):
    return dictionary_SNode_43[(SNode_42, GOAL_53, IDENTIFY39, SNode_43)]

dictionary_RESOLVE38 = {'true': 1.0, 'false': 1.0}

def f_RESOLVE38(RESOLVE38):
    return dictionary_RESOLVE38[RESOLVE38]

dictionary_GOAL_66 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_GOAL_66(SNode_43, GOAL_62, RESOLVE38, GOAL_66):
    return dictionary_GOAL_66[(SNode_43, GOAL_62, RESOLVE38, GOAL_66)]

dictionary_SNode_67 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_67(GOAL_66, NEED36, SNode_67):
    return dictionary_SNode_67[(GOAL_66, NEED36, SNode_67)]

dictionary_IDENTIFY41 = {'true': 1.0, 'false': 1.0}

def f_IDENTIFY41(IDENTIFY41):
    return dictionary_IDENTIFY41[IDENTIFY41]

dictionary_SNode_54 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_54(GOAL_53, IDENTIFY41, SNode_54):
    return dictionary_SNode_54[(GOAL_53, IDENTIFY41, SNode_54)]

dictionary_RESOLVE40 = {'true': 1.0, 'false': 1.0}

def f_RESOLVE40(RESOLVE40):
    return dictionary_RESOLVE40[RESOLVE40]

dictionary_GOAL_69 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_GOAL_69(SNode_54, GOAL_62, RESOLVE40, GOAL_69):
    return dictionary_GOAL_69[(SNode_54, GOAL_62, RESOLVE40, GOAL_69)]

dictionary_SNode_70 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_70(GOAL_69, NEED36, SNode_70):
    return dictionary_SNode_70[(GOAL_69, NEED36, SNode_70)]

dictionary_IDENTIFY43 = {'true': 1.0, 'false': 1.0}

def f_IDENTIFY43(IDENTIFY43):
    return dictionary_IDENTIFY43[IDENTIFY43]

dictionary_SNode_55 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_SNode_55(SNode_34, GOAL_53, IDENTIFY43, SNode_55):
    return dictionary_SNode_55[(SNode_34, GOAL_53, IDENTIFY43, SNode_55)]

dictionary_RESOLVE42 = {'true': 1.0, 'false': 1.0}

def f_RESOLVE42(RESOLVE42):
    return dictionary_RESOLVE42[RESOLVE42]

dictionary_GOAL_72 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_GOAL_72(SNode_55, GOAL_62, RESOLVE42, GOAL_72):
    return dictionary_GOAL_72[(SNode_55, GOAL_62, RESOLVE42, GOAL_72)]

dictionary_SNode_73 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_73(GOAL_72, NEED36, SNode_73):
    return dictionary_SNode_73[(GOAL_72, NEED36, SNode_73)]

dictionary_KINE29 = {'true': 1.0, 'false': 1.0}

def f_KINE29(KINE29):
    return dictionary_KINE29[KINE29]

dictionary_SNode_74 = {('false', 'true', 'false', 'true', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'true', 'true', 'true', 'true', 'false'): 1.0, ('true', 'true', 'true', 'true', 'true', 'false', 'true'): 1.0, ('false', 'false', 'false', 'false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true', 'false', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'false', 'false', 'true', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'false', 'false', 'true', 'false'): 1.0, ('false', 'false', 'true', 'false', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'true', 'false', 'false', 'true'): 1.0, ('false', 'false', 'true', 'true', 'true', 'false', 'false'): 1.0, ('false', 'true', 'false', 'true', 'true', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true', 'true', 'true', 'false'): 1.0, ('true', 'false', 'true', 'true', 'false', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'true', 'false', 'false'): 1.0, ('true', 'true', 'false', 'true', 'true', 'true', 'true'): 1.0, ('false', 'true', 'false', 'true', 'false', 'true', 'true'): 1.0, ('true', 'true', 'true', 'false', 'false', 'true', 'true'): 1.0, ('true', 'false', 'false', 'true', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'true', 'true', 'false'): 1.0, ('true', 'true', 'false', 'true', 'false', 'false', 'true'): 1.0, ('false', 'true', 'true', 'true', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'false', 'false', 'false'): 1.0, ('false', 'true', 'true', 'true', 'true', 'true', 'false'): 1.0, ('true', 'true', 'true', 'false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'false', 'false', 'true', 'false', 'false'): 1.0, ('false', 'false', 'false', 'true', 'true', 'false', 'false'): 1.0, ('true', 'true', 'true', 'false', 'false', 'true', 'false'): 1.0, ('false', 'false', 'true', 'false', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'false', 'false', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'false', 'true', 'false', 'false'): 1.0, ('true', 'true', 'false', 'false', 'true', 'false', 'true'): 1.0, ('false', 'false', 'true', 'true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'true', 'true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'true', 'true', 'true', 'false', 'true', 'true'): 1.0, ('true', 'true', 'false', 'true', 'false', 'true', 'true'): 1.0, ('true', 'false', 'true', 'true', 'true', 'false', 'true'): 1.0, ('false', 'false', 'false', 'true', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true', 'true', 'false', 'true'): 1.0, ('false', 'false', 'false', 'false', 'false', 'false', 'true'): 1.0, ('true', 'false', 'true', 'false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'true', 'false', 'false', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'false', 'false', 'true', 'true'): 1.0, ('true', 'true', 'true', 'true', 'false', 'true', 'false'): 1.0, ('false', 'false', 'false', 'true', 'false', 'true', 'false'): 1.0, ('true', 'false', 'true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'true', 'false', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'true', 'true', 'false', 'false', 'false'): 1.0, ('false', 'false', 'false', 'false', 'true', 'false', 'true'): 1.0, ('true', 'false', 'false', 'true', 'true', 'false', 'false'): 1.0, ('true', 'false', 'false', 'false', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'false', 'true', 'false', 'true', 'true'): 1.0, ('true', 'false', 'true', 'true', 'false', 'true', 'false'): 1.0, ('true', 'false', 'true', 'true', 'false', 'true', 'true'): 1.0, ('true', 'false', 'false', 'true', 'false', 'true', 'false'): 1.0, ('true', 'true', 'true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false', 'false', 'false'): 1.0, ('true', 'true', 'false', 'true', 'false', 'false', 'false'): 1.0, ('false', 'false', 'true', 'false', 'false', 'true', 'false'): 1.0, ('true', 'true', 'false', 'true', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false', 'true', 'false'): 1.0, ('false', 'false', 'true', 'true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'false', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false', 'false', 'true', 'true'): 1.0, ('true', 'false', 'false', 'true', 'false', 'false', 'true'): 1.0, ('true', 'false', 'false', 'false', 'true', 'true', 'true'): 1.0, ('false', 'true', 'true', 'false', 'false', 'true', 'true'): 1.0, ('true', 'false', 'true', 'false', 'true', 'true', 'true'): 1.0, ('true', 'false', 'true', 'false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false', 'false', 'true', 'true'): 1.0, ('true', 'true', 'false', 'false', 'true', 'false', 'false'): 1.0, ('true', 'true', 'true', 'true', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'true', 'false', 'false', 'false'): 1.0, ('true', 'true', 'true', 'true', 'false', 'false', 'true'): 1.0, ('true', 'false', 'false', 'true', 'true', 'true', 'false'): 1.0, ('true', 'true', 'false', 'true', 'true', 'false', 'false'): 1.0, ('false', 'true', 'false', 'false', 'false', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false', 'false', 'false', 'true'): 1.0, ('false', 'true', 'true', 'true', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false', 'true', 'false', 'true'): 1.0, ('true', 'false', 'true', 'false', 'false', 'false', 'true'): 1.0, ('false', 'true', 'true', 'true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'true', 'false', 'false', 'false'): 1.0, ('true', 'true', 'true', 'true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false', 'false', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false', 'false', 'true', 'false'): 1.0, ('true', 'true', 'true', 'false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true', 'true', 'true', 'true'): 1.0, ('true', 'true', 'false', 'false', 'true', 'true', 'true'): 1.0, ('true', 'false', 'true', 'false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'false', 'false', 'true'): 1.0, ('false', 'true', 'false', 'true', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true', 'false', 'true', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true', 'false', 'true'): 1.0, ('false', 'false', 'true', 'true', 'true', 'false', 'true'): 1.0, ('true', 'false', 'false', 'false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true', 'true', 'true', 'true'): 1.0, ('false', 'false', 'false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'true', 'true', 'true', 'false', 'true', 'false'): 1.0, ('true', 'false', 'true', 'false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'false', 'true', 'true', 'true'): 1.0, ('false', 'true', 'false', 'false', 'true', 'true', 'false'): 1.0, ('true', 'true', 'true', 'false', 'true', 'true', 'true'): 1.0, ('true', 'true', 'false', 'true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true', 'true', 'true', 'true'): 1.0, ('true', 'true', 'false', 'false', 'false', 'true', 'false'): 1.0, ('false', 'false', 'true', 'false', 'false', 'false', 'false'): 1.0}
def f_SNode_74(GOAL_62, SNode_64, SNode_67, SNode_70, SNode_73, KINE29, SNode_74):
    return dictionary_SNode_74[(GOAL_62, SNode_64, SNode_67, SNode_70, SNode_73, KINE29, SNode_74)]

dictionary_VECTOR44 = {'true': 1.0, 'false': 1.0}

def f_VECTOR44(VECTOR44):
    return dictionary_VECTOR44[VECTOR44]

dictionary_SNode_75 = {('true', 'true', 'false', 'true', 'true'): 1.0, ('true', 'true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true', 'true'): 1.0, ('true', 'true', 'false', 'false', 'true'): 1.0, ('false', 'true', 'true', 'true', 'true'): 1.0, ('false', 'true', 'true', 'true', 'false'): 1.0, ('false', 'true', 'false', 'false', 'true'): 1.0, ('false', 'false', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'false'): 1.0, ('false', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'false', 'true', 'false'): 1.0, ('true', 'true', 'true', 'true', 'true'): 1.0, ('false', 'true', 'false', 'true', 'false'): 1.0, ('true', 'false', 'true', 'false', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true'): 1.0, ('true', 'false', 'false', 'false', 'false'): 1.0, ('true', 'true', 'false', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'false', 'false', 'true'): 1.0, ('true', 'false', 'false', 'false', 'true'): 1.0, ('false', 'false', 'true', 'true', 'true'): 1.0, ('true', 'true', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true', 'true'): 1.0, ('true', 'true', 'true', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'false', 'false', 'false'): 1.0}
def f_SNode_75(SNode_4, GOAL_72, SNode_73, VECTOR44, SNode_75):
    return dictionary_SNode_75[(SNode_4, GOAL_72, SNode_73, VECTOR44, SNode_75)]

dictionary_EQUATION28 = {'true': 1.0, 'false': 1.0}

def f_EQUATION28(EQUATION28):
    return dictionary_EQUATION28[EQUATION28]

dictionary_GOAL_79 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_GOAL_79(SNode_74, SNode_75, EQUATION28, GOAL_79):
    return dictionary_GOAL_79[(SNode_74, SNode_75, EQUATION28, GOAL_79)]

dictionary_RApp5 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_RApp5(VECTOR27, GOAL_79, RApp5):
    return dictionary_RApp5[(VECTOR27, GOAL_79, RApp5)]

dictionary_GOAL_80 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_GOAL_80(SNode_75, EQUATION28, GOAL_80):
    return dictionary_GOAL_80[(SNode_75, EQUATION28, GOAL_80)]

dictionary_RApp6 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_RApp6(COMPO16, GOAL_80, RApp6):
    return dictionary_RApp6[(COMPO16, GOAL_80, RApp6)]

dictionary_GOAL_81 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_GOAL_81(RApp5, RApp6, GOAL_81):
    return dictionary_GOAL_81[(RApp5, RApp6, GOAL_81)]

dictionary_TRY25 = {'true': 1.0, 'false': 1.0}

def f_TRY25(TRY25):
    return dictionary_TRY25[TRY25]

dictionary_TRY24 = {('false', 'true'): 1.0, ('true', 'false'): 1.0, ('false', 'false'): 1.0, ('true', 'true'): 1.0}
def f_TRY24(TRY25, TRY24):
    return dictionary_TRY24[(TRY25, TRY24)]

dictionary_GOAL_83 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_GOAL_83(GOAL_81, TRY24, GOAL_83):
    return dictionary_GOAL_83[(GOAL_81, TRY24, GOAL_83)]

dictionary_CHOOSE47 = {'true': 1.0, 'false': 1.0}

def f_CHOOSE47(CHOOSE47):
    return dictionary_CHOOSE47[CHOOSE47]

dictionary_GOAL_84 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_GOAL_84(GOAL_83, CHOOSE47, GOAL_84):
    return dictionary_GOAL_84[(GOAL_83, CHOOSE47, GOAL_84)]

dictionary_SYSTEM46 = {'true': 1.0, 'false': 1.0}

def f_SYSTEM46(SYSTEM46):
    return dictionary_SYSTEM46[SYSTEM46]

dictionary_SNode_86 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_86(GOAL_84, SYSTEM46, SNode_86):
    return dictionary_SNode_86[(GOAL_84, SYSTEM46, SNode_86)]

dictionary_NEWTONS45 = {'true': 1.0, 'false': 1.0}

def f_NEWTONS45(NEWTONS45):
    return dictionary_NEWTONS45[NEWTONS45]

dictionary_SNode_156 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_156(SNode_86, NEWTONS45, SNode_156):
    return dictionary_SNode_156[(SNode_86, NEWTONS45, SNode_156)]

dictionary_DEFINE23 = {'true': 1.0, 'false': 1.0}

def f_DEFINE23(DEFINE23):
    return dictionary_DEFINE23[DEFINE23]

dictionary_GOAL_98 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_GOAL_98(GOAL_83, SNode_156, DEFINE23, GOAL_98):
    return dictionary_GOAL_98[(GOAL_83, SNode_156, DEFINE23, GOAL_98)]

dictionary_IDENTIFY22 = {'true': 1.0, 'false': 1.0}

def f_IDENTIFY22(IDENTIFY22):
    return dictionary_IDENTIFY22[IDENTIFY22]

dictionary_SNode_37 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_37(GOAL_98, IDENTIFY22, SNode_37):
    return dictionary_SNode_37[(GOAL_98, IDENTIFY22, SNode_37)]

dictionary_TRY26 = {('false', 'true'): 1.0, ('true', 'false'): 1.0, ('false', 'false'): 1.0, ('true', 'true'): 1.0}
def f_TRY26(TRY25, TRY26):
    return dictionary_TRY26[(TRY25, TRY26)]

dictionary_SNode_38 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_38(SNode_37, VAR20, SNode_38):
    return dictionary_SNode_38[(SNode_37, VAR20, SNode_38)]

dictionary_SNode_40 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_40(SNode_38, VALUE3, SNode_40):
    return dictionary_SNode_40[(SNode_38, VALUE3, SNode_40)]

dictionary_SNode_44 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_44(SNode_43, VAR20, SNode_44):
    return dictionary_SNode_44[(SNode_43, VAR20, SNode_44)]

dictionary_SNode_46 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_46(SNode_44, VALUE3, SNode_46):
    return dictionary_SNode_46[(SNode_44, VALUE3, SNode_46)]

dictionary_NULL48 = {'true': 1.0, 'false': 1.0}

def f_NULL48(NULL48):
    return dictionary_NULL48[NULL48]

dictionary_SNode_65 = {('true', 'true', 'false', 'true', 'true'): 1.0, ('true', 'true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true', 'true'): 1.0, ('true', 'true', 'false', 'false', 'true'): 1.0, ('false', 'true', 'true', 'true', 'true'): 1.0, ('false', 'true', 'true', 'true', 'false'): 1.0, ('false', 'true', 'false', 'false', 'true'): 1.0, ('false', 'false', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'false'): 1.0, ('false', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'false', 'true', 'false'): 1.0, ('true', 'true', 'true', 'true', 'true'): 1.0, ('false', 'true', 'false', 'true', 'false'): 1.0, ('true', 'false', 'true', 'false', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true'): 1.0, ('true', 'false', 'false', 'false', 'false'): 1.0, ('true', 'true', 'false', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'false', 'false', 'true'): 1.0, ('true', 'false', 'false', 'false', 'true'): 1.0, ('false', 'false', 'true', 'true', 'true'): 1.0, ('true', 'true', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true', 'true'): 1.0, ('true', 'true', 'true', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'false', 'false', 'false'): 1.0}
def f_SNode_65(SNode_29, GOAL_63, SNode_64, NULL48, SNode_65):
    return dictionary_SNode_65[(SNode_29, GOAL_63, SNode_64, NULL48, SNode_65)]

dictionary_SNode_68 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_SNode_68(GOAL_66, SNode_67, VECTOR44, SNode_68):
    return dictionary_SNode_68[(GOAL_66, SNode_67, VECTOR44, SNode_68)]

dictionary_SNode_71 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_SNode_71(GOAL_69, SNode_70, VECTOR44, SNode_71):
    return dictionary_SNode_71[(GOAL_69, SNode_70, VECTOR44, SNode_71)]

dictionary_FIND49 = {'true': 1.0, 'false': 1.0}

def f_FIND49(FIND49):
    return dictionary_FIND49[FIND49]

dictionary_GOAL_87 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_GOAL_87(GOAL_83, SNode_156, FIND49, GOAL_87):
    return dictionary_GOAL_87[(GOAL_83, SNode_156, FIND49, GOAL_87)]

dictionary_NORMAL50 = {'true': 1.0, 'false': 1.0}

def f_NORMAL50(NORMAL50):
    return dictionary_NORMAL50[NORMAL50]

dictionary_SNode_88 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_SNode_88(SNode_25, GOAL_87, NORMAL50, SNode_88):
    return dictionary_SNode_88[(SNode_25, GOAL_87, NORMAL50, SNode_88)]

dictionary_STRAT_90 = {'SNode_92_1': 1.0, 'SNode_91_2': 1.0}

def f_STRAT_90(STRAT_90):
    return dictionary_STRAT_90[STRAT_90]

dictionary_NORMAL52 = {'true': 1.0, 'false': 1.0}

def f_NORMAL52(NORMAL52):
    return dictionary_NORMAL52[NORMAL52]

dictionary_INCLINE51 = {('false', 'true'): 1.0, ('true', 'false'): 1.0, ('false', 'false'): 1.0, ('true', 'true'): 1.0}
def f_INCLINE51(NORMAL52, INCLINE51):
    return dictionary_INCLINE51[(NORMAL52, INCLINE51)]

dictionary_SNode_91 = {('false', 'false', 'false', 'SNode_91_2', 'false', 'true'): 1.0, ('false', 'false', 'false', 'SNode_91_2', 'true', 'false'): 1.0, ('false', 'true', 'true', 'SNode_91_2', 'true', 'false'): 1.0, ('false', 'false', 'false', 'SNode_91_2', 'false', 'false'): 1.0, ('false', 'false', 'true', 'SNode_92_1', 'false', 'false'): 1.0, ('true', 'false', 'false', 'SNode_92_1', 'false', 'false'): 1.0, ('false', 'true', 'false', 'SNode_91_2', 'true', 'true'): 1.0, ('false', 'false', 'true', 'SNode_91_2', 'false', 'false'): 1.0, ('true', 'false', 'false', 'SNode_91_2', 'true', 'false'): 1.0, ('false', 'true', 'true', 'SNode_92_1', 'false', 'true'): 1.0, ('true', 'true', 'false', 'SNode_92_1', 'false', 'false'): 1.0, ('false', 'true', 'true', 'SNode_91_2', 'true', 'true'): 1.0, ('true', 'true', 'false', 'SNode_91_2', 'true', 'false'): 1.0, ('true', 'false', 'false', 'SNode_91_2', 'false', 'true'): 1.0, ('false', 'false', 'true', 'SNode_92_1', 'false', 'true'): 1.0, ('true', 'true', 'true', 'SNode_92_1', 'true', 'false'): 1.0, ('false', 'true', 'false', 'SNode_91_2', 'true', 'false'): 1.0, ('true', 'false', 'true', 'SNode_91_2', 'false', 'true'): 1.0, ('true', 'false', 'false', 'SNode_91_2', 'true', 'true'): 1.0, ('false', 'true', 'true', 'SNode_92_1', 'true', 'true'): 1.0, ('true', 'true', 'false', 'SNode_92_1', 'true', 'true'): 1.0, ('true', 'true', 'true', 'SNode_91_2', 'true', 'true'): 1.0, ('true', 'false', 'true', 'SNode_91_2', 'false', 'false'): 1.0, ('true', 'true', 'false', 'SNode_92_1', 'false', 'true'): 1.0, ('true', 'false', 'true', 'SNode_92_1', 'true', 'true'): 1.0, ('true', 'false', 'false', 'SNode_91_2', 'false', 'false'): 1.0, ('false', 'false', 'false', 'SNode_91_2', 'true', 'true'): 1.0, ('false', 'true', 'false', 'SNode_92_1', 'false', 'true'): 1.0, ('true', 'true', 'true', 'SNode_91_2', 'false', 'true'): 1.0, ('true', 'false', 'false', 'SNode_92_1', 'false', 'true'): 1.0, ('true', 'false', 'true', 'SNode_92_1', 'false', 'false'): 1.0, ('false', 'true', 'true', 'SNode_92_1', 'true', 'false'): 1.0, ('true', 'true', 'false', 'SNode_92_1', 'true', 'false'): 1.0, ('true', 'true', 'false', 'SNode_91_2', 'false', 'true'): 1.0, ('true', 'false', 'true', 'SNode_91_2', 'true', 'false'): 1.0, ('true', 'true', 'false', 'SNode_91_2', 'false', 'false'): 1.0, ('false', 'true', 'false', 'SNode_92_1', 'false', 'false'): 1.0, ('false', 'false', 'true', 'SNode_92_1', 'true', 'true'): 1.0, ('true', 'false', 'true', 'SNode_91_2', 'true', 'true'): 1.0, ('true', 'false', 'true', 'SNode_92_1', 'false', 'true'): 1.0, ('true', 'true', 'true', 'SNode_91_2', 'false', 'false'): 1.0, ('false', 'false', 'true', 'SNode_91_2', 'true', 'false'): 1.0, ('true', 'true', 'true', 'SNode_92_1', 'false', 'true'): 1.0, ('false', 'true', 'false', 'SNode_91_2', 'false', 'true'): 1.0, ('true', 'true', 'true', 'SNode_92_1', 'true', 'true'): 1.0, ('false', 'true', 'false', 'SNode_92_1', 'true', 'true'): 1.0, ('false', 'false', 'false', 'SNode_92_1', 'true', 'false'): 1.0, ('false', 'true', 'false', 'SNode_92_1', 'true', 'false'): 1.0, ('false', 'true', 'true', 'SNode_91_2', 'false', 'true'): 1.0, ('false', 'false', 'false', 'SNode_92_1', 'true', 'true'): 1.0, ('true', 'true', 'true', 'SNode_91_2', 'true', 'false'): 1.0, ('false', 'true', 'true', 'SNode_91_2', 'false', 'false'): 1.0, ('true', 'true', 'false', 'SNode_91_2', 'true', 'true'): 1.0, ('false', 'false', 'false', 'SNode_92_1', 'false', 'true'): 1.0, ('true', 'false', 'false', 'SNode_92_1', 'true', 'false'): 1.0, ('true', 'false', 'true', 'SNode_92_1', 'true', 'false'): 1.0, ('false', 'true', 'false', 'SNode_91_2', 'false', 'false'): 1.0, ('false', 'false', 'true', 'SNode_92_1', 'true', 'false'): 1.0, ('false', 'true', 'true', 'SNode_92_1', 'false', 'false'): 1.0, ('true', 'true', 'true', 'SNode_92_1', 'false', 'false'): 1.0, ('false', 'false', 'true', 'SNode_91_2', 'true', 'true'): 1.0, ('true', 'false', 'false', 'SNode_92_1', 'true', 'true'): 1.0, ('false', 'false', 'true', 'SNode_91_2', 'false', 'true'): 1.0, ('false', 'false', 'false', 'SNode_92_1', 'false', 'false'): 1.0}
def f_SNode_91(SNode_88, SNode_12, SNode_13, STRAT_90, INCLINE51, SNode_91):
    return dictionary_SNode_91[(SNode_88, SNode_12, SNode_13, STRAT_90, INCLINE51, SNode_91)]

dictionary_HORIZ53 = {('false', 'true'): 1.0, ('true', 'false'): 1.0, ('false', 'false'): 1.0, ('true', 'true'): 1.0}
def f_HORIZ53(NORMAL52, HORIZ53):
    return dictionary_HORIZ53[(NORMAL52, HORIZ53)]

dictionary_BUGGY54 = {('false', 'true'): 1.0, ('true', 'false'): 1.0, ('false', 'false'): 1.0, ('true', 'true'): 1.0}
def f_BUGGY54(NORMAL52, BUGGY54):
    return dictionary_BUGGY54[(NORMAL52, BUGGY54)]

dictionary_SNode_92 = {('true', 'SNode_92_1', 'false', 'false'): 1.0, ('true', 'SNode_91_2', 'false', 'true'): 1.0, ('true', 'SNode_92_1', 'false', 'true'): 1.0, ('true', 'SNode_92_1', 'true', 'true'): 1.0, ('true', 'SNode_91_2', 'true', 'true'): 1.0, ('false', 'SNode_92_1', 'true', 'false'): 1.0, ('false', 'SNode_91_2', 'false', 'true'): 1.0, ('true', 'SNode_92_1', 'true', 'false'): 1.0, ('false', 'SNode_92_1', 'false', 'true'): 1.0, ('true', 'SNode_91_2', 'true', 'false'): 1.0, ('true', 'SNode_91_2', 'false', 'false'): 1.0, ('false', 'SNode_91_2', 'false', 'false'): 1.0, ('false', 'SNode_92_1', 'false', 'false'): 1.0, ('false', 'SNode_91_2', 'true', 'true'): 1.0, ('false', 'SNode_92_1', 'true', 'true'): 1.0, ('false', 'SNode_91_2', 'true', 'false'): 1.0}
def f_SNode_92(SNode_12, STRAT_90, BUGGY54, SNode_92):
    return dictionary_SNode_92[(SNode_12, STRAT_90, BUGGY54, SNode_92)]

dictionary_IDENTIFY55 = {'true': 1.0, 'false': 1.0}

def f_IDENTIFY55(IDENTIFY55):
    return dictionary_IDENTIFY55[IDENTIFY55]

dictionary_SNode_93 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_SNode_93(GOAL_87, SNode_88, IDENTIFY55, SNode_93):
    return dictionary_SNode_93[(GOAL_87, SNode_88, IDENTIFY55, SNode_93)]

dictionary_WEIGHT56 = {'true': 1.0, 'false': 1.0}

def f_WEIGHT56(WEIGHT56):
    return dictionary_WEIGHT56[WEIGHT56]

dictionary_SNode_94 = {('true', 'true', 'false', 'true', 'true'): 1.0, ('true', 'true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true', 'true'): 1.0, ('true', 'true', 'false', 'false', 'true'): 1.0, ('false', 'true', 'true', 'true', 'true'): 1.0, ('false', 'true', 'true', 'true', 'false'): 1.0, ('false', 'true', 'false', 'false', 'true'): 1.0, ('false', 'false', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'false'): 1.0, ('false', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'false', 'true', 'false'): 1.0, ('true', 'true', 'true', 'true', 'true'): 1.0, ('false', 'true', 'false', 'true', 'false'): 1.0, ('true', 'false', 'true', 'false', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true'): 1.0, ('true', 'false', 'false', 'false', 'false'): 1.0, ('true', 'true', 'false', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'false', 'false', 'true'): 1.0, ('true', 'false', 'false', 'false', 'true'): 1.0, ('false', 'false', 'true', 'true', 'true'): 1.0, ('true', 'true', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true', 'true'): 1.0, ('true', 'true', 'true', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'false', 'false', 'false'): 1.0}
def f_SNode_94(SNode_16, SNode_33, GOAL_87, WEIGHT56, SNode_94):
    return dictionary_SNode_94[(SNode_16, SNode_33, GOAL_87, WEIGHT56, SNode_94)]

dictionary_WEIGHT57 = {'true': 1.0, 'false': 1.0}

def f_WEIGHT57(WEIGHT57):
    return dictionary_WEIGHT57[WEIGHT57]

dictionary_SNode_95 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_95(SNode_94, WEIGHT57, SNode_95):
    return dictionary_SNode_95[(SNode_94, WEIGHT57, SNode_95)]

dictionary_SNode_97 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_SNode_97(GOAL_87, SNode_94, IDENTIFY55, SNode_97):
    return dictionary_SNode_97[(GOAL_87, SNode_94, IDENTIFY55, SNode_97)]

dictionary_FIND58 = {'true': 1.0, 'false': 1.0}

def f_FIND58(FIND58):
    return dictionary_FIND58[FIND58]

dictionary_GOAL_99 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_GOAL_99(GOAL_98, FIND58, GOAL_99):
    return dictionary_GOAL_99[(GOAL_98, FIND58, GOAL_99)]

dictionary_IDENTIFY59 = {'true': 1.0, 'false': 1.0}

def f_IDENTIFY59(IDENTIFY59):
    return dictionary_IDENTIFY59[IDENTIFY59]

dictionary_SNode_100 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_100(GOAL_98, IDENTIFY59, SNode_100):
    return dictionary_SNode_100[(GOAL_98, IDENTIFY59, SNode_100)]

dictionary_FORCE60 = {'true': 1.0, 'false': 1.0}

def f_FORCE60(FORCE60):
    return dictionary_FORCE60[FORCE60]

dictionary_SNode_102 = {('true', 'true', 'false', 'true', 'true'): 1.0, ('true', 'true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true', 'true'): 1.0, ('true', 'true', 'false', 'false', 'true'): 1.0, ('false', 'true', 'true', 'true', 'true'): 1.0, ('false', 'true', 'true', 'true', 'false'): 1.0, ('false', 'true', 'false', 'false', 'true'): 1.0, ('false', 'false', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'false'): 1.0, ('false', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'false', 'true', 'false'): 1.0, ('true', 'true', 'true', 'true', 'true'): 1.0, ('false', 'true', 'false', 'true', 'false'): 1.0, ('true', 'false', 'true', 'false', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true'): 1.0, ('true', 'false', 'false', 'false', 'false'): 1.0, ('true', 'true', 'false', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'false', 'false', 'true'): 1.0, ('true', 'false', 'false', 'false', 'true'): 1.0, ('false', 'false', 'true', 'true', 'true'): 1.0, ('true', 'true', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true', 'true'): 1.0, ('true', 'true', 'true', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'false', 'false', 'false'): 1.0}
def f_SNode_102(GOAL_87, SNode_88, SNode_94, FORCE60, SNode_102):
    return dictionary_SNode_102[(GOAL_87, SNode_88, SNode_94, FORCE60, SNode_102)]

dictionary_APPLY61 = {'true': 1.0, 'false': 1.0}

def f_APPLY61(APPLY61):
    return dictionary_APPLY61[APPLY61]

dictionary_GOAL_103 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_GOAL_103(GOAL_83, SNode_102, APPLY61, GOAL_103):
    return dictionary_GOAL_103[(GOAL_83, SNode_102, APPLY61, GOAL_103)]

dictionary_CHOOSE62 = {'true': 1.0, 'false': 1.0}

def f_CHOOSE62(CHOOSE62):
    return dictionary_CHOOSE62[CHOOSE62]

dictionary_GOAL_104 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_GOAL_104(GOAL_103, CHOOSE62, GOAL_104):
    return dictionary_GOAL_104[(GOAL_103, CHOOSE62, GOAL_104)]

dictionary_SNode_106 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_106(GOAL_104, MAXIMIZE34, SNode_106):
    return dictionary_SNode_106[(GOAL_104, MAXIMIZE34, SNode_106)]

dictionary_SNode_152 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_152(SNode_106, AXIS33, SNode_152):
    return dictionary_SNode_152[(SNode_106, AXIS33, SNode_152)]

dictionary_WRITE63 = {'true': 1.0, 'false': 1.0}

def f_WRITE63(WRITE63):
    return dictionary_WRITE63[WRITE63]

dictionary_GOAL_107 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_GOAL_107(GOAL_103, SNode_152, WRITE63, GOAL_107):
    return dictionary_GOAL_107[(GOAL_103, SNode_152, WRITE63, GOAL_107)]

dictionary_WRITE64 = {'true': 1.0, 'false': 1.0}

def f_WRITE64(WRITE64):
    return dictionary_WRITE64[WRITE64]

dictionary_GOAL_108 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_GOAL_108(GOAL_107, WRITE64, GOAL_108):
    return dictionary_GOAL_108[(GOAL_107, WRITE64, GOAL_108)]

dictionary_GOAL_109 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_GOAL_109(GOAL_107, WRITE64, GOAL_109):
    return dictionary_GOAL_109[(GOAL_107, WRITE64, GOAL_109)]

dictionary_GOAL65 = {'true': 1.0, 'false': 1.0}

def f_GOAL65(GOAL65):
    return dictionary_GOAL65[GOAL65]

dictionary_GOAL_110 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_GOAL_110(GOAL_109, GOAL65, GOAL_110):
    return dictionary_GOAL_110[(GOAL_109, GOAL65, GOAL_110)]

dictionary_GOAL66 = {'true': 1.0, 'false': 1.0}

def f_GOAL66(GOAL66):
    return dictionary_GOAL66[GOAL66]

dictionary_GOAL_111 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_GOAL_111(GOAL_109, GOAL66, GOAL_111):
    return dictionary_GOAL_111[(GOAL_109, GOAL66, GOAL_111)]

dictionary_NEED67 = {'true': 1.0, 'false': 1.0}

def f_NEED67(NEED67):
    return dictionary_NEED67[NEED67]

dictionary_RApp7 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_RApp7(NEED67, GOAL_109, RApp7):
    return dictionary_RApp7[(NEED67, GOAL_109, RApp7)]

dictionary_RApp8 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_RApp8(NEED36, GOAL_111, RApp8):
    return dictionary_RApp8[(NEED36, GOAL_111, RApp8)]

dictionary_SNode_112 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_112(RApp7, RApp8, SNode_112):
    return dictionary_SNode_112[(RApp7, RApp8, SNode_112)]

dictionary_GOAL68 = {'true': 1.0, 'false': 1.0}

def f_GOAL68(GOAL68):
    return dictionary_GOAL68[GOAL68]

dictionary_GOAL_113 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_GOAL_113(GOAL_110, GOAL68, GOAL_113):
    return dictionary_GOAL_113[(GOAL_110, GOAL68, GOAL_113)]

dictionary_GOAL_114 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_GOAL_114(GOAL_110, GOAL68, GOAL_114):
    return dictionary_GOAL_114[(GOAL_110, GOAL68, GOAL_114)]

dictionary_SNode_115 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_115(GOAL_114, NEED36, SNode_115):
    return dictionary_SNode_115[(GOAL_114, NEED36, SNode_115)]

dictionary_VECTOR69 = {'true': 1.0, 'false': 1.0}

def f_VECTOR69(VECTOR69):
    return dictionary_VECTOR69[VECTOR69]

dictionary_SNode_116 = {('true', 'true', 'true', 'true', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'false', 'true', 'true'): 1.0, ('false', 'true', 'false', 'false', 'false', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true', 'true'): 1.0, ('false', 'true', 'false', 'true', 'false', 'true'): 1.0, ('false', 'true', 'true', 'true', 'true', 'false'): 1.0, ('true', 'true', 'true', 'true', 'true', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'false', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'false', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true', 'true', 'true'): 1.0, ('false', 'false', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'false', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'true', 'true', 'true', 'true'): 1.0, ('true', 'true', 'true', 'false', 'true', 'true'): 1.0, ('true', 'false', 'true', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'true', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true', 'false', 'false'): 1.0, ('false', 'false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true', 'true', 'true'): 1.0, ('true', 'true', 'true', 'false', 'true', 'false'): 1.0, ('true', 'false', 'true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false', 'false', 'false'): 1.0, ('false', 'false', 'false', 'false', 'true', 'false'): 1.0, ('true', 'true', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false', 'false', 'true'): 1.0, ('false', 'false', 'false', 'true', 'true', 'false'): 1.0, ('true', 'true', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false', 'true'): 1.0, ('true', 'false', 'true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'true', 'true', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'false', 'true', 'true'): 1.0, ('true', 'true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'true', 'false', 'false', 'false', 'true'): 1.0, ('true', 'false', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true', 'false', 'true'): 1.0, ('false', 'false', 'false', 'true', 'false', 'false'): 1.0, ('true', 'true', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false', 'false'): 1.0, ('false', 'false', 'true', 'false', 'true', 'false'): 1.0, ('false', 'false', 'false', 'true', 'true', 'true'): 1.0, ('true', 'false', 'true', 'true', 'true', 'true'): 1.0, ('true', 'true', 'false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false', 'false', 'true'): 1.0, ('true', 'false', 'false', 'true', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'true', 'false'): 1.0}
def f_SNode_116(SNode_95, SNode_97, GOAL_114, SNode_115, VECTOR69, SNode_116):
    return dictionary_SNode_116[(SNode_95, SNode_97, GOAL_114, SNode_115, VECTOR69, SNode_116)]

dictionary_SNode_117 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_117(GOAL_113, NEED36, SNode_117):
    return dictionary_SNode_117[(GOAL_113, NEED36, SNode_117)]

dictionary_VECTOR70 = {'true': 1.0, 'false': 1.0}

def f_VECTOR70(VECTOR70):
    return dictionary_VECTOR70[VECTOR70]

dictionary_SNode_118 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_SNode_118(SNode_91, GOAL_113, VECTOR70, SNode_118):
    return dictionary_SNode_118[(SNode_91, GOAL_113, VECTOR70, SNode_118)]

dictionary_EQUAL71 = {'true': 1.0, 'false': 1.0}

def f_EQUAL71(EQUAL71):
    return dictionary_EQUAL71[EQUAL71]

dictionary_SNode_119 = {('true', 'true', 'false', 'true', 'true'): 1.0, ('true', 'true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true', 'true'): 1.0, ('true', 'true', 'false', 'false', 'true'): 1.0, ('false', 'true', 'true', 'true', 'true'): 1.0, ('false', 'true', 'true', 'true', 'false'): 1.0, ('false', 'true', 'false', 'false', 'true'): 1.0, ('false', 'false', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'false'): 1.0, ('false', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'false', 'true', 'false'): 1.0, ('true', 'true', 'true', 'true', 'true'): 1.0, ('false', 'true', 'false', 'true', 'false'): 1.0, ('true', 'false', 'true', 'false', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true'): 1.0, ('true', 'false', 'false', 'false', 'false'): 1.0, ('true', 'true', 'false', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'false', 'false', 'true'): 1.0, ('true', 'false', 'false', 'false', 'true'): 1.0, ('false', 'false', 'true', 'true', 'true'): 1.0, ('true', 'true', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true', 'true'): 1.0, ('true', 'true', 'true', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'false', 'false', 'false'): 1.0}
def f_SNode_119(SNode_93, SNode_117, SNode_118, EQUAL71, SNode_119):
    return dictionary_SNode_119[(SNode_93, SNode_117, SNode_118, EQUAL71, SNode_119)]

dictionary_SNode_120 = {('true', 'true', 'true', 'true', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'false', 'true', 'true'): 1.0, ('false', 'true', 'false', 'false', 'false', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true', 'true'): 1.0, ('false', 'true', 'false', 'true', 'false', 'true'): 1.0, ('false', 'true', 'true', 'true', 'true', 'false'): 1.0, ('true', 'true', 'true', 'true', 'true', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'false', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'false', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true', 'true', 'true'): 1.0, ('false', 'false', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'false', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'true', 'true', 'true', 'true'): 1.0, ('true', 'true', 'true', 'false', 'true', 'true'): 1.0, ('true', 'false', 'true', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'true', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true', 'false', 'false'): 1.0, ('false', 'false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true', 'true', 'true'): 1.0, ('true', 'true', 'true', 'false', 'true', 'false'): 1.0, ('true', 'false', 'true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false', 'false', 'false'): 1.0, ('false', 'false', 'false', 'false', 'true', 'false'): 1.0, ('true', 'true', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false', 'false', 'true'): 1.0, ('false', 'false', 'false', 'true', 'true', 'false'): 1.0, ('true', 'true', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false', 'true'): 1.0, ('true', 'false', 'true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'true', 'true', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'false', 'true', 'true'): 1.0, ('true', 'true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'true', 'false', 'false', 'false', 'true'): 1.0, ('true', 'false', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true', 'false', 'true'): 1.0, ('false', 'false', 'false', 'true', 'false', 'false'): 1.0, ('true', 'true', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false', 'false'): 1.0, ('false', 'false', 'true', 'false', 'true', 'false'): 1.0, ('false', 'false', 'false', 'true', 'true', 'true'): 1.0, ('true', 'false', 'true', 'true', 'true', 'true'): 1.0, ('true', 'true', 'false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false', 'false', 'true'): 1.0, ('true', 'false', 'false', 'true', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'true', 'false'): 1.0}
def f_SNode_120(SNode_92, SNode_93, GOAL_113, SNode_117, VECTOR69, SNode_120):
    return dictionary_SNode_120[(SNode_92, SNode_93, GOAL_113, SNode_117, VECTOR69, SNode_120)]

dictionary_GOAL72 = {'true': 1.0, 'false': 1.0}

def f_GOAL72(GOAL72):
    return dictionary_GOAL72[GOAL72]

dictionary_GOAL_121 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_GOAL_121(GOAL_109, GOAL72, GOAL_121):
    return dictionary_GOAL_121[(GOAL_109, GOAL72, GOAL_121)]

dictionary_SNode_122 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_122(GOAL_121, NEED36, SNode_122):
    return dictionary_SNode_122[(GOAL_121, NEED36, SNode_122)]

dictionary_VECTOR73 = {'true': 1.0, 'false': 1.0}

def f_VECTOR73(VECTOR73):
    return dictionary_VECTOR73[VECTOR73]

dictionary_SNode_123 = {('true', 'true', 'true', 'true', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'false', 'true', 'true'): 1.0, ('false', 'true', 'false', 'false', 'false', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true', 'true'): 1.0, ('false', 'true', 'false', 'true', 'false', 'true'): 1.0, ('false', 'true', 'true', 'true', 'true', 'false'): 1.0, ('true', 'true', 'true', 'true', 'true', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'false', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'false', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true', 'true', 'true'): 1.0, ('false', 'false', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'false', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'true', 'true', 'true', 'true'): 1.0, ('true', 'true', 'true', 'false', 'true', 'true'): 1.0, ('true', 'false', 'true', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'true', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true', 'false', 'false'): 1.0, ('false', 'false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true', 'true', 'true'): 1.0, ('true', 'true', 'true', 'false', 'true', 'false'): 1.0, ('true', 'false', 'true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false', 'false', 'false'): 1.0, ('false', 'false', 'false', 'false', 'true', 'false'): 1.0, ('true', 'true', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false', 'false', 'true'): 1.0, ('false', 'false', 'false', 'true', 'true', 'false'): 1.0, ('true', 'true', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false', 'true'): 1.0, ('true', 'false', 'true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'true', 'true', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'false', 'true', 'true'): 1.0, ('true', 'true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'true', 'false', 'false', 'false', 'true'): 1.0, ('true', 'false', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true', 'false', 'true'): 1.0, ('false', 'false', 'false', 'true', 'false', 'false'): 1.0, ('true', 'true', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false', 'false'): 1.0, ('false', 'false', 'true', 'false', 'true', 'false'): 1.0, ('false', 'false', 'false', 'true', 'true', 'true'): 1.0, ('true', 'false', 'true', 'true', 'true', 'true'): 1.0, ('true', 'true', 'false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false', 'false', 'true'): 1.0, ('true', 'false', 'false', 'true', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'true', 'false'): 1.0}
def f_SNode_123(SNode_4, SNode_100, GOAL_121, SNode_122, VECTOR73, SNode_123):
    return dictionary_SNode_123[(SNode_4, SNode_100, GOAL_121, SNode_122, VECTOR73, SNode_123)]

dictionary_NEWTONS74 = {'true': 1.0, 'false': 1.0}

def f_NEWTONS74(NEWTONS74):
    return dictionary_NEWTONS74[NEWTONS74]

dictionary_SNode_124 = {('true', 'true', 'true', 'true', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'false', 'true', 'true'): 1.0, ('false', 'true', 'false', 'false', 'false', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true', 'true'): 1.0, ('false', 'true', 'false', 'true', 'false', 'true'): 1.0, ('false', 'true', 'true', 'true', 'true', 'false'): 1.0, ('true', 'true', 'true', 'true', 'true', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'false', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'false', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true', 'true', 'true'): 1.0, ('false', 'false', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'false', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'true', 'true', 'true', 'true'): 1.0, ('true', 'true', 'true', 'false', 'true', 'true'): 1.0, ('true', 'false', 'true', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'true', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true', 'false', 'false'): 1.0, ('false', 'false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true', 'true', 'true'): 1.0, ('true', 'true', 'true', 'false', 'true', 'false'): 1.0, ('true', 'false', 'true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false', 'false', 'false'): 1.0, ('false', 'false', 'false', 'false', 'true', 'false'): 1.0, ('true', 'true', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false', 'false', 'true'): 1.0, ('false', 'false', 'false', 'true', 'true', 'false'): 1.0, ('true', 'true', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false', 'true'): 1.0, ('true', 'false', 'true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'true', 'true', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'false', 'true', 'true'): 1.0, ('true', 'true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'true', 'false', 'false', 'false', 'true'): 1.0, ('true', 'false', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true', 'false', 'true'): 1.0, ('false', 'false', 'false', 'true', 'false', 'false'): 1.0, ('true', 'true', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false', 'false'): 1.0, ('false', 'false', 'true', 'false', 'true', 'false'): 1.0, ('false', 'false', 'false', 'true', 'true', 'true'): 1.0, ('true', 'false', 'true', 'true', 'true', 'true'): 1.0, ('true', 'true', 'false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false', 'false', 'true'): 1.0, ('true', 'false', 'false', 'true', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'true', 'false'): 1.0}
def f_SNode_124(SNode_37, GOAL_109, SNode_112, SNode_122, NEWTONS74, SNode_124):
    return dictionary_SNode_124[(SNode_37, GOAL_109, SNode_112, SNode_122, NEWTONS74, SNode_124)]

dictionary_SUM75 = {'true': 1.0, 'false': 1.0}

def f_SUM75(SUM75):
    return dictionary_SUM75[SUM75]

dictionary_SNode_125 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_SNode_125(GOAL_109, SNode_112, SUM75, SNode_125):
    return dictionary_SNode_125[(GOAL_109, SNode_112, SUM75, SNode_125)]

dictionary_GOAL_126 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_GOAL_126(GOAL_108, GOAL65, GOAL_126):
    return dictionary_GOAL_126[(GOAL_108, GOAL65, GOAL_126)]

dictionary_GOAL_127 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_GOAL_127(GOAL_108, GOAL66, GOAL_127):
    return dictionary_GOAL_127[(GOAL_108, GOAL66, GOAL_127)]

dictionary_RApp9 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_RApp9(NEED67, GOAL_108, RApp9):
    return dictionary_RApp9[(NEED67, GOAL_108, RApp9)]

dictionary_RApp10 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_RApp10(NEED36, GOAL_127, RApp10):
    return dictionary_RApp10[(NEED36, GOAL_127, RApp10)]

dictionary_SNode_128 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_128(RApp9, RApp10, SNode_128):
    return dictionary_SNode_128[(RApp9, RApp10, SNode_128)]

dictionary_GOAL_129 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_GOAL_129(GOAL_126, GOAL68, GOAL_129):
    return dictionary_GOAL_129[(GOAL_126, GOAL68, GOAL_129)]

dictionary_GOAL_130 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_GOAL_130(GOAL_126, GOAL68, GOAL_130):
    return dictionary_GOAL_130[(GOAL_126, GOAL68, GOAL_130)]

dictionary_SNode_131 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_131(GOAL_130, NEED36, SNode_131):
    return dictionary_SNode_131[(GOAL_130, NEED36, SNode_131)]

dictionary_SNode_132 = {('true', 'true', 'true', 'true', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'false', 'true', 'true'): 1.0, ('false', 'true', 'false', 'false', 'false', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true', 'true'): 1.0, ('false', 'true', 'false', 'true', 'false', 'true'): 1.0, ('false', 'true', 'true', 'true', 'true', 'false'): 1.0, ('true', 'true', 'true', 'true', 'true', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'false', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'false', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true', 'true', 'true'): 1.0, ('false', 'false', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'false', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'true', 'true', 'true', 'true'): 1.0, ('true', 'true', 'true', 'false', 'true', 'true'): 1.0, ('true', 'false', 'true', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'true', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true', 'false', 'false'): 1.0, ('false', 'false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true', 'true', 'true'): 1.0, ('true', 'true', 'true', 'false', 'true', 'false'): 1.0, ('true', 'false', 'true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false', 'false', 'false'): 1.0, ('false', 'false', 'false', 'false', 'true', 'false'): 1.0, ('true', 'true', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false', 'false', 'true'): 1.0, ('false', 'false', 'false', 'true', 'true', 'false'): 1.0, ('true', 'true', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false', 'true'): 1.0, ('true', 'false', 'true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'true', 'true', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'false', 'true', 'true'): 1.0, ('true', 'true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'true', 'false', 'false', 'false', 'true'): 1.0, ('true', 'false', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true', 'false', 'true'): 1.0, ('false', 'false', 'false', 'true', 'false', 'false'): 1.0, ('true', 'true', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false', 'false'): 1.0, ('false', 'false', 'true', 'false', 'true', 'false'): 1.0, ('false', 'false', 'false', 'true', 'true', 'true'): 1.0, ('true', 'false', 'true', 'true', 'true', 'true'): 1.0, ('true', 'true', 'false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false', 'false', 'true'): 1.0, ('true', 'false', 'false', 'true', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'true', 'false'): 1.0}
def f_SNode_132(SNode_95, SNode_97, GOAL_130, SNode_131, VECTOR69, SNode_132):
    return dictionary_SNode_132[(SNode_95, SNode_97, GOAL_130, SNode_131, VECTOR69, SNode_132)]

dictionary_SNode_133 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_133(GOAL_129, NEED36, SNode_133):
    return dictionary_SNode_133[(GOAL_129, NEED36, SNode_133)]

dictionary_SNode_134 = {('true', 'true', 'true', 'true', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'false', 'true', 'true'): 1.0, ('false', 'true', 'false', 'false', 'false', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true', 'true'): 1.0, ('false', 'true', 'false', 'true', 'false', 'true'): 1.0, ('false', 'true', 'true', 'true', 'true', 'false'): 1.0, ('true', 'true', 'true', 'true', 'true', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'false', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'false', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true', 'true', 'true'): 1.0, ('false', 'false', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'false', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'true', 'true', 'true', 'true'): 1.0, ('true', 'true', 'true', 'false', 'true', 'true'): 1.0, ('true', 'false', 'true', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'true', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true', 'false', 'false'): 1.0, ('false', 'false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true', 'true', 'true'): 1.0, ('true', 'true', 'true', 'false', 'true', 'false'): 1.0, ('true', 'false', 'true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false', 'false', 'false'): 1.0, ('false', 'false', 'false', 'false', 'true', 'false'): 1.0, ('true', 'true', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false', 'false', 'true'): 1.0, ('false', 'false', 'false', 'true', 'true', 'false'): 1.0, ('true', 'true', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false', 'true'): 1.0, ('true', 'false', 'true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'true', 'true', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'false', 'true', 'true'): 1.0, ('true', 'true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'true', 'false', 'false', 'false', 'true'): 1.0, ('true', 'false', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true', 'false', 'true'): 1.0, ('false', 'false', 'false', 'true', 'false', 'false'): 1.0, ('true', 'true', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false', 'false'): 1.0, ('false', 'false', 'true', 'false', 'true', 'false'): 1.0, ('false', 'false', 'false', 'true', 'true', 'true'): 1.0, ('true', 'false', 'true', 'true', 'true', 'true'): 1.0, ('true', 'true', 'false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false', 'false', 'true'): 1.0, ('true', 'false', 'false', 'true', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'true', 'false'): 1.0}
def f_SNode_134(SNode_91, SNode_93, GOAL_129, SNode_133, VECTOR73, SNode_134):
    return dictionary_SNode_134[(SNode_91, SNode_93, GOAL_129, SNode_133, VECTOR73, SNode_134)]

dictionary_SNode_135 = {('true', 'true', 'true', 'true', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'false', 'true', 'true'): 1.0, ('false', 'true', 'false', 'false', 'false', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true', 'true'): 1.0, ('false', 'true', 'false', 'true', 'false', 'true'): 1.0, ('false', 'true', 'true', 'true', 'true', 'false'): 1.0, ('true', 'true', 'true', 'true', 'true', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'false', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'false', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true', 'true', 'true'): 1.0, ('false', 'false', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'false', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'true', 'true', 'true', 'true'): 1.0, ('true', 'true', 'true', 'false', 'true', 'true'): 1.0, ('true', 'false', 'true', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'true', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true', 'false', 'false'): 1.0, ('false', 'false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true', 'true', 'true'): 1.0, ('true', 'true', 'true', 'false', 'true', 'false'): 1.0, ('true', 'false', 'true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false', 'false', 'false'): 1.0, ('false', 'false', 'false', 'false', 'true', 'false'): 1.0, ('true', 'true', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false', 'false', 'true'): 1.0, ('false', 'false', 'false', 'true', 'true', 'false'): 1.0, ('true', 'true', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false', 'true'): 1.0, ('true', 'false', 'true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'true', 'true', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'false', 'true', 'true'): 1.0, ('true', 'true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'true', 'false', 'false', 'false', 'true'): 1.0, ('true', 'false', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true', 'false', 'true'): 1.0, ('false', 'false', 'false', 'true', 'false', 'false'): 1.0, ('true', 'true', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false', 'false'): 1.0, ('false', 'false', 'true', 'false', 'true', 'false'): 1.0, ('false', 'false', 'false', 'true', 'true', 'true'): 1.0, ('true', 'false', 'true', 'true', 'true', 'true'): 1.0, ('true', 'true', 'false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false', 'false', 'true'): 1.0, ('true', 'false', 'false', 'true', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'true', 'false'): 1.0}
def f_SNode_135(SNode_92, SNode_93, GOAL_129, SNode_133, VECTOR69, SNode_135):
    return dictionary_SNode_135[(SNode_92, SNode_93, GOAL_129, SNode_133, VECTOR69, SNode_135)]

dictionary_SNode_154 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_154(GOAL_121, NEED36, SNode_154):
    return dictionary_SNode_154[(GOAL_121, NEED36, SNode_154)]

dictionary_SNode_136 = {('true', 'true', 'true', 'true', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'false', 'true', 'true'): 1.0, ('false', 'true', 'false', 'false', 'false', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true', 'true'): 1.0, ('false', 'true', 'false', 'true', 'false', 'true'): 1.0, ('false', 'true', 'true', 'true', 'true', 'false'): 1.0, ('true', 'true', 'true', 'true', 'true', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'false', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'false', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true', 'true', 'true'): 1.0, ('false', 'false', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'false', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'true', 'true', 'true', 'true'): 1.0, ('true', 'true', 'true', 'false', 'true', 'true'): 1.0, ('true', 'false', 'true', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'true', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true', 'false', 'false'): 1.0, ('false', 'false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true', 'true', 'true'): 1.0, ('true', 'true', 'true', 'false', 'true', 'false'): 1.0, ('true', 'false', 'true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false', 'false', 'false'): 1.0, ('false', 'false', 'false', 'false', 'true', 'false'): 1.0, ('true', 'true', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false', 'false', 'true'): 1.0, ('false', 'false', 'false', 'true', 'true', 'false'): 1.0, ('true', 'true', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false', 'true'): 1.0, ('true', 'false', 'true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'true', 'true', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'false', 'true', 'true'): 1.0, ('true', 'true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'true', 'false', 'false', 'false', 'true'): 1.0, ('true', 'false', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true', 'false', 'true'): 1.0, ('false', 'false', 'false', 'true', 'false', 'false'): 1.0, ('true', 'true', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false', 'false'): 1.0, ('false', 'false', 'true', 'false', 'true', 'false'): 1.0, ('false', 'false', 'false', 'true', 'true', 'true'): 1.0, ('true', 'false', 'true', 'true', 'true', 'true'): 1.0, ('true', 'true', 'false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false', 'false', 'true'): 1.0, ('true', 'false', 'false', 'true', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'true', 'false'): 1.0}
def f_SNode_136(SNode_37, GOAL_108, SNode_128, SNode_154, NEWTONS74, SNode_136):
    return dictionary_SNode_136[(SNode_37, GOAL_108, SNode_128, SNode_154, NEWTONS74, SNode_136)]

dictionary_SNode_137 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_SNode_137(GOAL_108, SNode_128, SUM75, SNode_137):
    return dictionary_SNode_137[(GOAL_108, SNode_128, SUM75, SNode_137)]

dictionary_GOAL_142 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_GOAL_142(SNode_116, SNode_125, EQUATION28, GOAL_142):
    return dictionary_GOAL_142[(SNode_116, SNode_125, EQUATION28, GOAL_142)]

dictionary_GOAL_143 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_GOAL_143(SNode_116, SNode_132, EQUATION28, GOAL_143):
    return dictionary_GOAL_143[(SNode_116, SNode_132, EQUATION28, GOAL_143)]

dictionary_GOAL_146 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_GOAL_146(SNode_132, SNode_137, EQUATION28, GOAL_146):
    return dictionary_GOAL_146[(SNode_132, SNode_137, EQUATION28, GOAL_146)]

dictionary_RApp11 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_RApp11(VECTOR27, GOAL_142, RApp11):
    return dictionary_RApp11[(VECTOR27, GOAL_142, RApp11)]

dictionary_RApp12 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_RApp12(COMPO16, GOAL_143, RApp12):
    return dictionary_RApp12[(COMPO16, GOAL_143, RApp12)]

dictionary_RApp13 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_RApp13(VECTOR27, GOAL_146, RApp13):
    return dictionary_RApp13[(VECTOR27, GOAL_146, RApp13)]

dictionary_GOAL_147 = {('true', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true'): 1.0}
def f_GOAL_147(RApp11, RApp12, RApp13, GOAL_147):
    return dictionary_GOAL_147[(RApp11, RApp12, RApp13, GOAL_147)]

dictionary_TRY76 = {'true': 1.0, 'false': 1.0}

def f_TRY76(TRY76):
    return dictionary_TRY76[TRY76]

dictionary_GOAL_149 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_GOAL_149(GOAL_147, TRY76, GOAL_149):
    return dictionary_GOAL_149[(GOAL_147, TRY76, GOAL_149)]

dictionary_APPLY77 = {'true': 1.0, 'false': 1.0}

def f_APPLY77(APPLY77):
    return dictionary_APPLY77[APPLY77]

dictionary_GOAL_150 = {('true', 'true', 'false', 'true', 'true'): 1.0, ('true', 'true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true', 'true'): 1.0, ('true', 'true', 'false', 'false', 'true'): 1.0, ('false', 'true', 'true', 'true', 'true'): 1.0, ('false', 'true', 'true', 'true', 'false'): 1.0, ('false', 'true', 'false', 'false', 'true'): 1.0, ('false', 'false', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'false'): 1.0, ('false', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'false', 'true', 'false'): 1.0, ('true', 'true', 'true', 'true', 'true'): 1.0, ('false', 'true', 'false', 'true', 'false'): 1.0, ('true', 'false', 'true', 'false', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true'): 1.0, ('true', 'false', 'false', 'false', 'false'): 1.0, ('true', 'true', 'false', 'true', 'false'): 1.0, ('false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false'): 1.0, ('false', 'false', 'false', 'false', 'true'): 1.0, ('true', 'false', 'false', 'false', 'true'): 1.0, ('false', 'false', 'true', 'true', 'true'): 1.0, ('true', 'true', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true', 'true'): 1.0, ('true', 'true', 'true', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true', 'false'): 1.0, ('false', 'true', 'false', 'false', 'false'): 1.0}
def f_GOAL_150(SNode_20, SNode_37, GOAL_149, APPLY77, GOAL_150):
    return dictionary_GOAL_150[(SNode_20, SNode_37, GOAL_149, APPLY77, GOAL_150)]

dictionary_GRAV78 = {'true': 1.0, 'false': 1.0}

def f_GRAV78(GRAV78):
    return dictionary_GRAV78[GRAV78]

dictionary_SNode_151 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_SNode_151(GOAL_150, GRAV78, SNode_151):
    return dictionary_SNode_151[(GOAL_150, GRAV78, SNode_151)]

dictionary_GOAL_153 = {('true', 'true', 'false'): 1.0, ('false', 'false', 'true'): 1.0, ('true', 'true', 'true'): 1.0, ('true', 'false', 'false'): 1.0, ('false', 'false', 'false'): 1.0, ('false', 'true', 'false'): 1.0, ('false', 'true', 'true'): 1.0, ('true', 'false', 'true'): 1.0}
def f_GOAL_153(GOAL_108, GOAL72, GOAL_153):
    return dictionary_GOAL_153[(GOAL_108, GOAL72, GOAL_153)]

dictionary_SNode_155 = {('true', 'true', 'true', 'true', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false', 'true', 'false'): 1.0, ('true', 'false', 'false', 'false', 'true', 'true'): 1.0, ('false', 'true', 'false', 'false', 'false', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true', 'true'): 1.0, ('false', 'true', 'false', 'true', 'false', 'true'): 1.0, ('false', 'true', 'true', 'true', 'true', 'false'): 1.0, ('true', 'true', 'true', 'true', 'true', 'true'): 1.0, ('false', 'true', 'true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'false', 'false', 'false', 'false'): 1.0, ('false', 'false', 'true', 'true', 'false', 'true'): 1.0, ('true', 'true', 'false', 'true', 'true', 'true'): 1.0, ('true', 'false', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'false', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true', 'true', 'true'): 1.0, ('false', 'false', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'false', 'false', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'true', 'true', 'true', 'false'): 1.0, ('false', 'false', 'true', 'true', 'true', 'true'): 1.0, ('true', 'true', 'true', 'false', 'true', 'true'): 1.0, ('true', 'false', 'true', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'true', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'true', 'true', 'false'): 1.0, ('false', 'true', 'false', 'true', 'false', 'false'): 1.0, ('false', 'false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'true', 'true', 'true'): 1.0, ('true', 'true', 'true', 'false', 'true', 'false'): 1.0, ('true', 'false', 'true', 'false', 'false', 'true'): 1.0, ('true', 'true', 'false', 'false', 'false', 'false'): 1.0, ('false', 'false', 'false', 'false', 'true', 'false'): 1.0, ('true', 'true', 'false', 'true', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false', 'true', 'true'): 1.0, ('true', 'false', 'false', 'false', 'false', 'true'): 1.0, ('false', 'false', 'false', 'true', 'true', 'false'): 1.0, ('true', 'true', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'false', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false', 'true'): 1.0, ('true', 'false', 'true', 'false', 'false', 'false'): 1.0, ('true', 'false', 'true', 'false', 'true', 'false'): 1.0, ('false', 'true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'true', 'true', 'false', 'false', 'false'): 1.0, ('false', 'true', 'false', 'false', 'true', 'true'): 1.0, ('true', 'true', 'true', 'true', 'false', 'false'): 1.0, ('true', 'true', 'false', 'false', 'false', 'true'): 1.0, ('true', 'false', 'false', 'true', 'true', 'false'): 1.0, ('true', 'false', 'true', 'true', 'false', 'true'): 1.0, ('false', 'true', 'false', 'true', 'true', 'true'): 1.0, ('false', 'false', 'false', 'true', 'false', 'true'): 1.0, ('false', 'false', 'false', 'true', 'false', 'false'): 1.0, ('true', 'true', 'false', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'false', 'false'): 1.0, ('false', 'true', 'true', 'false', 'false', 'false'): 1.0, ('false', 'false', 'true', 'false', 'true', 'false'): 1.0, ('false', 'false', 'false', 'true', 'true', 'true'): 1.0, ('true', 'false', 'true', 'true', 'true', 'true'): 1.0, ('true', 'true', 'false', 'true', 'false', 'true'): 1.0, ('true', 'true', 'true', 'false', 'false', 'true'): 1.0, ('true', 'false', 'false', 'true', 'false', 'true'): 1.0, ('false', 'false', 'true', 'false', 'true', 'true'): 1.0, ('false', 'false', 'true', 'true', 'true', 'false'): 1.0}
def f_SNode_155(SNode_4, SNode_100, GOAL_153, SNode_154, VECTOR44, SNode_155):
    return dictionary_SNode_155[(SNode_4, SNode_100, GOAL_153, SNode_154, VECTOR44, SNode_155)]

functions = [f_GOAL_2, f_SNode_3, f_SNode_4, f_SNode_5, f_SNode_6, f_SNode_7, f_DISPLACEM0, f_RApp1, f_GIVEN_1, f_RApp2, f_SNode_8, f_SNode_9, f_SNode_10, f_SNode_11, f_SNode_12, f_SNode_13, f_SNode_14, f_SNode_15, f_SNode_16, f_SNode_17, f_SNode_18, f_SNode_19, f_NEED1, f_SNode_20, f_GRAV2, f_SNode_21, f_VALUE3, f_SNode_24, f_SLIDING4, f_SNode_25, f_CONSTANT5, f_SNode_26, f_KNOWN6, f_VELOCITY7, f_SNode_47, f_RApp3, f_KNOWN8, f_RApp4, f_SNode_27, f_COMPO16, f_GOAL_48, f_TRY12, f_TRY11, f_GOAL_49, f_CHOOSE19, f_GOAL_50, f_SYSTEM18, f_SNode_51, f_KINEMATI17, f_SNode_52, f_IDENTIFY10, f_GOAL_53, f_IDENTIFY9, f_SNode_28, f_TRY13, f_TRY14, f_TRY15, f_VAR20, f_SNode_29, f_SNode_31, f_GIVEN21, f_SNode_33, f_SNode_34, f_VECTOR27, f_APPLY32, f_GOAL_56, f_CHOOSE35, f_GOAL_57, f_MAXIMIZE34, f_SNode_59, f_AXIS33, f_SNode_60, f_WRITE31, f_GOAL_61, f_WRITE30, f_GOAL_62, f_RESOLVE37, f_GOAL_63, f_NEED36, f_SNode_64, f_SNode_41, f_SNode_42, f_IDENTIFY39, f_SNode_43, f_RESOLVE38, f_GOAL_66, f_SNode_67, f_IDENTIFY41, f_SNode_54, f_RESOLVE40, f_GOAL_69, f_SNode_70, f_IDENTIFY43, f_SNode_55, f_RESOLVE42, f_GOAL_72, f_SNode_73, f_KINE29, f_SNode_74, f_VECTOR44, f_SNode_75, f_EQUATION28, f_GOAL_79, f_RApp5, f_GOAL_80, f_RApp6, f_GOAL_81, f_TRY25, f_TRY24, f_GOAL_83, f_CHOOSE47, f_GOAL_84, f_SYSTEM46, f_SNode_86, f_NEWTONS45, f_SNode_156, f_DEFINE23, f_GOAL_98, f_IDENTIFY22, f_SNode_37, f_TRY26, f_SNode_38, f_SNode_40, f_SNode_44, f_SNode_46, f_NULL48, f_SNode_65, f_SNode_68, f_SNode_71, f_FIND49, f_GOAL_87, f_NORMAL50, f_SNode_88, f_STRAT_90, f_NORMAL52, f_INCLINE51, f_SNode_91, f_HORIZ53, f_BUGGY54, f_SNode_92, f_IDENTIFY55, f_SNode_93, f_WEIGHT56, f_SNode_94, f_WEIGHT57, f_SNode_95, f_SNode_97, f_FIND58, f_GOAL_99, f_IDENTIFY59, f_SNode_100, f_FORCE60, f_SNode_102, f_APPLY61, f_GOAL_103, f_CHOOSE62, f_GOAL_104, f_SNode_106, f_SNode_152, f_WRITE63, f_GOAL_107, f_WRITE64, f_GOAL_108, f_GOAL_109, f_GOAL65, f_GOAL_110, f_GOAL66, f_GOAL_111, f_NEED67, f_RApp7, f_RApp8, f_SNode_112, f_GOAL68, f_GOAL_113, f_GOAL_114, f_SNode_115, f_VECTOR69, f_SNode_116, f_SNode_117, f_VECTOR70, f_SNode_118, f_EQUAL71, f_SNode_119, f_SNode_120, f_GOAL72, f_GOAL_121, f_SNode_122, f_VECTOR73, f_SNode_123, f_NEWTONS74, f_SNode_124, f_SUM75, f_SNode_125, f_GOAL_126, f_GOAL_127, f_RApp9, f_RApp10, f_SNode_128, f_GOAL_129, f_GOAL_130, f_SNode_131, f_SNode_132, f_SNode_133, f_SNode_134, f_SNode_135, f_SNode_154, f_SNode_136, f_SNode_137, f_GOAL_142, f_GOAL_143, f_GOAL_146, f_RApp11, f_RApp12, f_RApp13, f_GOAL_147, f_TRY76, f_GOAL_149, f_APPLY77, f_GOAL_150, f_GRAV78, f_SNode_151, f_GOAL_153, f_SNode_155]
domains_dict = {'GOAL68': ['false', 'true'], 'RESOLVE37': ['false', 'true'], 'FIND58': ['false', 'true'], 'NULL48': ['false', 'true'], 'SNode_137': ['false', 'true'], 'SNode_136': ['false', 'true'], 'SNode_135': ['false', 'true'], 'SNode_134': ['false', 'true'], 'RESOLVE38': ['false', 'true'], 'SNode_132': ['false', 'true'], 'SNode_131': ['false', 'true'], 'GOAL66': ['false', 'true'], 'WRITE31': ['false', 'true'], 'WRITE30': ['false', 'true'], 'VALUE3': ['false', 'true'], 'HORIZ53': ['false', 'true'], 'SNode_155': ['false', 'true'], 'NEED36': ['false', 'true'], 'SNode_33': ['false', 'true'], 'SNode_31': ['false', 'true'], 'SNode_37': ['false', 'true'], 'SNode_34': ['false', 'true'], 'SNode_38': ['false', 'true'], 'IDENTIFY59': ['false', 'true'], 'SNode_112': ['false', 'true'], 'VECTOR69': ['false', 'true'], 'SNode_120': ['false', 'true'], 'FIND49': ['false', 'true'], 'SNode_123': ['false', 'true'], 'SNode_124': ['false', 'true'], 'SNode_125': ['false', 'true'], 'SYSTEM46': ['false', 'true'], 'SNode_128': ['false', 'true'], 'GOAL65': ['false', 'true'], 'RESOLVE40': ['false', 'true'], 'RESOLVE42': ['false', 'true'], 'IDENTIFY10': ['false', 'true'], 'GOAL_150': ['false', 'true'], 'VAR20': ['false', 'true'], 'SNode_133': ['false', 'true'], 'INCLINE51': ['false', 'true'], 'IDENTIFY9': ['false', 'true'], 'COMPO16': ['false', 'true'], 'VELOCITY7': ['false', 'true'], 'GOAL_147': ['false', 'true'], 'SLIDING4': ['false', 'true'], 'SNode_88': ['false', 'true'], 'VECTOR27': ['false', 'true'], 'APPLY77': ['false', 'true'], 'DISPLACEM0': ['false', 'true'], 'KINE29': ['false', 'true'], 'GOAL_149': ['false', 'true'], 'SNode_154': ['false', 'true'], 'SNode_156': ['false', 'true'], 'SNode_151': ['false', 'true'], 'SNode_152': ['false', 'true'], 'IDENTIFY22': ['false', 'true'], 'GOAL_142': ['false', 'true'], 'EQUATION28': ['false', 'true'], 'GOAL_146': ['false', 'true'], 'SNode_19': ['false', 'true'], 'SNode_18': ['false', 'true'], 'VECTOR73': ['false', 'true'], 'SNode_11': ['false', 'true'], 'SNode_10': ['false', 'true'], 'SNode_13': ['false', 'true'], 'SNode_12': ['false', 'true'], 'SNode_15': ['false', 'true'], 'SNode_14': ['false', 'true'], 'SNode_17': ['false', 'true'], 'SNode_16': ['false', 'true'], 'SNode_74': ['false', 'true'], 'GIVEN_1': ['false', 'true'], 'SNode_91': ['false', 'true'], 'GOAL72': ['false', 'true'], 'SNode_93': ['false', 'true'], 'SNode_92': ['false', 'true'], 'SNode_95': ['false', 'true'], 'SNode_94': ['false', 'true'], 'SNode_97': ['false', 'true'], 'GOAL_49': ['false', 'true'], 'GOAL_48': ['false', 'true'], 'VECTOR70': ['false', 'true'], 'SYSTEM18': ['false', 'true'], 'GIVEN21': ['false', 'true'], 'SNode_60': ['false', 'true'], 'IDENTIFY39': ['false', 'true'], 'KNOWN6': ['false', 'true'], 'FORCE60': ['false', 'true'], 'KNOWN8': ['false', 'true'], 'CHOOSE62': ['false', 'true'], 'GOAL_130': ['false', 'true'], 'RApp8': ['false', 'true'], 'SNode_86': ['false', 'true'], 'SNode_68': ['false', 'true'], 'SNode_64': ['false', 'true'], 'RApp1': ['false', 'true'], 'RApp2': ['false', 'true'], 'RApp3': ['false', 'true'], 'RApp4': ['false', 'true'], 'RApp5': ['false', 'true'], 'RApp6': ['false', 'true'], 'RApp7': ['false', 'true'], 'RApp9': ['false', 'true'], 'DEFINE23': ['false', 'true'], 'GOAL_56': ['false', 'true'], 'GOAL_57': ['false', 'true'], 'GOAL_50': ['false', 'true'], 'GOAL_53': ['false', 'true'], 'CHOOSE47': ['false', 'true'], 'SNode_119': ['false', 'true'], 'MAXIMIZE34': ['false', 'true'], 'GOAL_129': ['false', 'true'], 'SNode_118': ['false', 'true'], 'GOAL_127': ['false', 'true'], 'GOAL_126': ['false', 'true'], 'KINEMATI17': ['false', 'true'], 'GOAL_121': ['false', 'true'], 'SNode_65': ['false', 'true'], 'IDENTIFY43': ['false', 'true'], 'SNode_75': ['false', 'true'], 'IDENTIFY41': ['false', 'true'], 'SNode_73': ['false', 'true'], 'SNode_71': ['false', 'true'], 'SNode_70': ['false', 'true'], 'NEED1': ['false', 'true'], 'WEIGHT56': ['false', 'true'], 'SNode_67': ['false', 'true'], 'TRY15': ['false', 'true'], 'TRY14': ['false', 'true'], 'TRY11': ['false', 'true'], 'TRY13': ['false', 'true'], 'TRY12': ['false', 'true'], 'GOAL_61': ['false', 'true'], 'GOAL_63': ['false', 'true'], 'GOAL_62': ['false', 'true'], 'GOAL_66': ['false', 'true'], 'GOAL_69': ['false', 'true'], 'SNode_115': ['false', 'true'], 'GOAL_143': ['false', 'true'], 'GOAL_113': ['false', 'true'], 'GOAL_110': ['false', 'true'], 'GOAL_111': ['false', 'true'], 'SNode_9': ['false', 'true'], 'SNode_8': ['false', 'true'], 'GOAL_114': ['false', 'true'], 'SNode_5': ['false', 'true'], 'CHOOSE19': ['false', 'true'], 'SNode_7': ['false', 'true'], 'SNode_6': ['false', 'true'], 'SNode_3': ['false', 'true'], 'SNode_42': ['false', 'true'], 'SNode_43': ['false', 'true'], 'SNode_40': ['false', 'true'], 'SNode_41': ['false', 'true'], 'SNode_46': ['false', 'true'], 'SNode_47': ['false', 'true'], 'SNode_44': ['false', 'true'], 'NORMAL52': ['false', 'true'], 'NORMAL50': ['false', 'true'], 'IDENTIFY55': ['false', 'true'], 'APPLY32': ['false', 'true'], 'GOAL_83': ['false', 'true'], 'SNode_54': ['false', 'true'], 'WRITE64': ['false', 'true'], 'WRITE63': ['false', 'true'], 'SUM75': ['false', 'true'], 'SNode_117': ['false', 'true'], 'GOAL_79': ['false', 'true'], 'GOAL_87': ['false', 'true'], 'GOAL_72': ['false', 'true'], 'GOAL_2': ['false', 'true'], 'GOAL_104': ['false', 'true'], 'GOAL_107': ['false', 'true'], 'NEWTONS45': ['false', 'true'], 'GOAL_103': ['false', 'true'], 'GOAL_84': ['false', 'true'], 'GOAL_109': ['false', 'true'], 'GOAL_108': ['false', 'true'], 'BUGGY54': ['false', 'true'], 'SNode_116': ['false', 'true'], 'SNode_55': ['false', 'true'], 'CONSTANT5': ['false', 'true'], 'GOAL_81': ['false', 'true'], 'GOAL_80': ['false', 'true'], 'SNode_51': ['false', 'true'], 'SNode_52': ['false', 'true'], 'SNode_59': ['false', 'true'], 'WEIGHT57': ['false', 'true'], 'RApp12': ['false', 'true'], 'RApp13': ['false', 'true'], 'RApp10': ['false', 'true'], 'RApp11': ['false', 'true'], 'GRAV78': ['false', 'true'], 'STRAT_90': ['SNode_92_1', 'SNode_91_2'], 'APPLY61': ['false', 'true'], 'SNode_4': ['false', 'true'], 'CHOOSE35': ['false', 'true'], 'SNode_102': ['false', 'true'], 'GRAV2': ['false', 'true'], 'SNode_100': ['false', 'true'], 'SNode_106': ['false', 'true'], 'TRY24': ['false', 'true'], 'TRY25': ['false', 'true'], 'TRY26': ['false', 'true'], 'AXIS33': ['false', 'true'], 'NEED67': ['false', 'true'], 'GOAL_98': ['false', 'true'], 'GOAL_99': ['false', 'true'], 'TRY76': ['false', 'true'], 'SNode_20': ['false', 'true'], 'SNode_21': ['false', 'true'], 'SNode_24': ['false', 'true'], 'SNode_25': ['false', 'true'], 'SNode_26': ['false', 'true'], 'SNode_27': ['false', 'true'], 'SNode_28': ['false', 'true'], 'SNode_29': ['false', 'true'], 'SNode_122': ['false', 'true'], 'EQUAL71': ['false', 'true'], 'NEWTONS74': ['false', 'true'], 'VECTOR44': ['false', 'true'], 'GOAL_153': ['false', 'true']}

def create_graph():
    g = build_graph(
        *functions,
        domains = domains_dict)
    g.name = 'A3-data/r'
    return g

def create_bbn():
    g = build_bbn(
        *functions,
        domains = domains_dict)
    g.name = 'A3-data/r'
    return g

