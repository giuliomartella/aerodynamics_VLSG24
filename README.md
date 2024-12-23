# aerodynamics_VLSG24
Aerodynamics course final project.
Contributors:
- Viola Monti
- Lorenzo Bascucci
- Simone Nannetti
- Giulio Martella

# Hess-Smith Method
We have completed the given `main` function.

## Theodorsen's Angle
1. Use `CommandGenXfoil_cp` to generate a `.txt` file with a sequence of Xfoil commands.
2. Place the generated file in Xfoil to obtain the distribution of `cp` as the angle of attack varies between -1.3° and -1°.
3. Use `data_manipulator` within the `alfa_cp` folder to create a `.mat` file.
4. Run `Theodorsen_cp` to calculate Theodorsen's angle using the `cp` distribution.
5. Run `Theodorsen_hs` to calculate Theodorsen's angle using the Hess-Smith method.
6. Run `MeanCamberLine` or `tryMean` to calculate Theodorsen's angle by defining the mean line in two different ways, then numerically integrating.

## Separation and Transition
1. Use `CommandGenXfoil` to generate a `.txt` document with Xfoil commands.
2. Place the generated file in Xfoil to obtain the `cf` behavior.

## Weissinger Method
Run `reportScript` within the `reportAnalysis` folder.
