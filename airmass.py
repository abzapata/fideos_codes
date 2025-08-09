from numpy import cos, pi, round

def airmass(x):
    AM = 1/(cos((90.-x)*pi/180.))
    return AM

ang = input("Enter altitud from the horizont in degress: ")  # Angle of elevation from the horizont
ang = float(ang)

res = round(airmass(ang),2)

print('Airmass = ', res)