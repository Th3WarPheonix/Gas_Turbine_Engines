
def newton1(fcn, dfcn, x1=0, precision=1e-6, **kwargs):
    """Classic newton raphson method; quadratic convergence"""
    x2 = x1 + precision*10
    while abs(x2-x1)>precision:
        x1 = x2
        x2 = x1-fcn(x1, **kwargs)/dfcn(x1, **kwargs)
    return x2

def newton2(fcn, x1=0, precision=1e-6, goal=0, **kwargs):
    """Does not require a derivative; quadratic convergence. If function-to-be-zeroed requires additional arguments, they will be passed in the order they are inputted"""
    delta = precision/10
    x2 = x1 + precision*10
    while abs(x2-x1)>precision:
        x1 = x2
        y1 = (fcn(x1+delta, **kwargs)-fcn(x1, **kwargs))/delta
        dx = (goal-fcn(x1, **kwargs))/y1
        x2 += dx
    return x2

def secant(fcn, xn1=0, precision=1e-6, **kwargs):
    """Secant method convergence rate of 1.618 (golden ratio); newtons method with finite difference"""
    xn2 = xn1 + precision*10
    while abs(xn2-xn1) > precision:
        xn3 = (xn2*fcn(xn1, **kwargs)-xn1*fcn(xn2, **kwargs))/(fcn(xn1, **kwargs)-fcn(xn2, **kwargs))
        xn1, xn2 = xn2, xn3
    return xn3

if __name__ == '__main__':
    initx = 1e5
    prec = 1e-10
    x1 = newton1(fcn1, dfcn1, x1=initx, precision=prec)
    x2 = newton2(fcn1, x1=initx, precision=prec)
    x3 = secant(fcn1, xn1=initx, precision=prec)
    print(x1, x2, x3)


