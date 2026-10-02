# COMPREHENSIVE NONLINEAR SEISMIC ASSESSMENT OF A POST-FIRE MULTI-DEGREE-FREEDOM WITH FOUR ELEMENTS STRUCTURE : AN OPENSEES FRAMEWORK FOR STATIC PUSHOVER, CYCLIC DEGRADATION, STATIC TIME HISTORY AND DYNAMIC TIME-HISTORY ANALYSIS

# ارزیابی جامع غیرخطی سازه چند درجه آزادی پسا آتش سوزی ستون فولادی : چارچوبی مبتنی بر اوپنسیس در پایتون برای تحلیل پوش‌آور استاتیکی، تحلیل چرخه‌ای همراه با کاهش سختی، تحلیل تاریخچه زمانی استاتیکی و دینامیکی، با در نظر گرفتن شاخص‌های آسیب شکل‌پذیری در المان و سازه و ارزیابی شاخص ظرفیت انرژی اتلاف‌شده


This Python Scripts encode the temperature-dependent degradation of structural steel per the experimental data.
It takes the steel temperature T (°C) and returns the reduction factors for elastic modulus (kE), yield strength (kFy), ultimate strength (kFu), and the associated ultimate strain (esu).
The stepwise logic reflects the distinct stages of material deterioration: a gradual decline up to 400 °C, a precipitous drop in stiffness and strength in the 400–500 °C range, and severe residual capacity beyond 700 °C.
These factors are critical for nonlinear thermo-mechanical analysis, enabling accurate prediction of load-bearing capacity, deflections, and failure modes during fire exposure.
Implementation follows the tabulated values in the referenced paper, providing a simple yet experimentally grounded fire-response material model for structural steel. 

![alt text](https://github.com/salardelavar/OPENSEES_SALAR/blob/main/EIGHT_ANALYSIS_DUCTILITY_DAMAGE_INDEX_%26_ENERGY_DISSIPATION_CAPACITY_INDEX/57_POST-FIRE/EXAMPLE_04_MDOF/COVER.png)

# EQUIVALENT SDOF SYSTEM DERIVATION VIA DISPLACEMENT-BASED SEISMIC DESIGN PROCEDURE WITH PUSHOVER AND FREE-VIBRATION ANALYSIS
![alt text](https://github.com/salardelavar/OPENSEES_SALAR/blob/main/EIGHT_ANALYSIS_DUCTILITY_DAMAGE_INDEX_%26_ENERGY_DISSIPATION_CAPACITY_INDEX/57_POST-FIRE/EXAMPLE_04_MDOF/COVER_DISPLACEMENT_BASED_PUSHOVER_%26_FREE_VIBRATION.png)
