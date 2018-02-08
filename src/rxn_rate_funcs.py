#Reaction rate functions

BT=5e-6 ##TODO

##TODO
def rate_forward(cca,ccl,ccacam):
    return cca*(BT-ccacam)

##TODO
def rate_backward(cca,ccl,ccacam):
    return ccacam