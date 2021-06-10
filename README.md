# Temperate-Mosquito-DDE

The code presented here forms the basis for the West Nile virus (WNV) model code which can be found at https://github.com/davewi13/WNV_model.  I would recommend in any future analyses that the WNV model code be used, as it is the most up-to-date.  By setting the transmission rates to zero this code can easily be used to model only mosquito populations if virus dynamics are not required.

For code used in the paper "Modelling the Effects of Temperature on Temperature Mosquito Seasonal Abundance" please see file Temperate_mosquito_DDE_code.f90.  However, please note that small modifications after publication mean that the code in the file "Chapter 2 DDE code.f90" should be used for all future analyses (results will be similar but not identical between the two sets of code).  The thesis which "Chapter 2 DDE code.f90" relates to can be found at http://theses.gla.ac.uk/8450/1/2017ewingphd.pdf.  Alternatively use the WNV model code with disease transmission turned off as specified above. 

For code used in the paper "Uncovering mechanisms behind mosquito seasonality by integrating mathematical models and daily empirical population data: Culex pipiens in the UK" please see the file Ewing et al 2018 - Data paper.f90.

For code used in my PhD please see the file named by the relevant chapter number.

The DDE solver code can be found at http://www.radford.edu/~thompson/ffddes/
