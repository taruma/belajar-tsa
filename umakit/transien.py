# Persamaan yang digunakan pada kuliah Transient

#if __name__ == '__main__':
#    import basicmath as bm
#else:
#    from . import basicmath as bm
    
from scipy import constants
g = constants.g
    
def koef_tunnel(lamda, panjang_tunnel, diameter_pipa, g=g):
    return lamda * panjang_tunnel / (2 * g * diameter_pipa)


def koef_total(koef_tanki, 
                    luas_tunnel, luas_tank, koef_tunnel):
    return (koef_tanki * pow(luas_tunnel / luas_tank, 2) 
                + koef_tunnel)


def head_loss(koefisien_tunnel, kecepatan_awal):
    return koefisien_tunnel * pow(kecepatan_awal, 2)