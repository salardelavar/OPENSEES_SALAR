# This document presents a comprehensive structural analysis of a reinforced concrete structure using OpenSees and Python which are written by Salar Delavar Ghashghaei (Qashqai). It includes various advanced procedures such as pushover, thermal, seismic, blast, and fire analyses, as well as optimization and uncertainty assessments. Please note that the content may not be entirely free of errors or inaccuracies.

Including:

[1] Pushover Analysis

[2] Hysteretic Pushover Analysis

[3] Structural Analysis Considering Sectional Strengthening with Steel Plates or FRP Composites

[4] Thermal Analysis Due to Fire Exposure

[5] Seismic Response Analysis Under Earthquake Loading

[6] Incremental Dynamic Analysis (IDA)

[7] Computation of Acceleration, Velocity, and Displacement Response Spectra

[8] Free Vibration (Modal) Analysis

[9] Structural Analysis Under Blast, Harmonic, and Wind-Induced Impact Loads

[10] Post-Buckling Behavior Analysis of Frame Columns

[11] Sensitivity Analysis of Reinforcement Ratio and Section Dimensions on Structural Ductility Ratio and Structural Behavior Coefficient

[12] Optimization of Rebar Diameter and Section Depth Based on Structural Ductility Ratio and Structural Behavior Coefficient and Structural Period

[13] Elastic Structural Analysis

[14] Structural Analysis Under Incrementally Increasing Distributed Loads

[15] Sequential Fire-Then-Earthquake Dynamic Analysis

[16] Assessment of the Structural Ductility Damage Index

[17] Uncertainty Analysis of RC Frames and Probabilistic Seismic Assessment

[18] Evaluation of Structural Ductility Damage Index Considering Pinned Beam Connections

[19] Sequential Explosion impact load-Then-Fire Analysis

[20] Analysis of Creep and Shrinkage Effects in Concrete Frames

[21] Soil–Structure Interaction with Foundation Consideration

[22] Nonlinear Static and Dynamic Analysis of Ultra-High Strength Concrete Frame Structure

[23] Progressive Collapse of a Concrete Frame Structure

[24] Utilizing Parallel Processing Procedures for the Simultaneous Execution of Nonlinear Static and Dynamic Concrete Structural Analysis

[25] Utilizing Parallel Processing Procedures for the Simultaneous Execution of Optimization Analysis

[26] Optimization of the concrete confinement coefficient in various column and beam sections through nonlinear static analysis, aiming to maximize the ductility ratio, using the Newton–Raphson algorithm. Also Sensitivity Analysis with Column Depth and Confiment Enhancement Ratio.

[27] Nonlinear Static and Dynamic Analysis of a Reinforced Concrete Frame with a Viscous Damper Using OpenSees and Python.

[28] Comparative Study of Elastic and Inelastic Structural Behavior through Pushover Dynamic Analysis.

[29] Comparative Study of Elastic and Inelastic Structural Behavior through Pushover Dynamic Analysis Due to Impact Load.

[30] Implementation of Response Spectrum Anlaysis with CQC Modal Combination.



تحلیل جامع یک سازه بتنی در نرم‌افزار اوپنسیس و پایتون شامل:

1- تحلیل پوش‌آور

2- پوش‌آور هیسترزی

3- تحلیل سازه با تاثیر تقویت مقطع بتنی با ورق فلزی یا اف آر پی

4- تحلیل حرارتی ناشی از آتش سوزی

5- تحلیل پاسخ لرزه‌ای ناشی از زلزله

6- تحلیل دینامیکی افزایشی

7- محاسبه طیف پاسخ شتاب، سرعت و جابجایی سازه

8- تحلیل ارتعاش آزاد سازه

9- تحلیل سازه با اثر بار ضربه ناشی از انفجار و هارمونیک و باد

10- تحلیل سازه با رفتار پس‌کمانشی ستون های قاب

11- تحلیل حساسیت میگرد و ابعاد سازه بر نسبت شکل پذیری و ضریب رفتار سازه

12- بهینه‌یابی قطر میلگرد و عمق مقطع براساس نسبت شکل پذیری و ضریب رفتار سازه و زمان تناوب سازه

13- تحلیل الاستیک سازه

14- تحلیل سازه با بارگسترده افزایشی

15- تحلیل حرارتی ناشی از آتش‌سوزی و سپس تحلیل دینامیکی ناشی از زلزله

16-محاسبه شاخص آسیب شکل‌پذیری سازه

17- تحلیل احتمالاتی آسیب لرزه ای سازه در شرایط عدم قطعیت

18- تحلیل شاخص آسیب شکل پذیری سازه با در نظر گرفتن اتصالات مفصلی تیر ها

19- تحلیل دینامیکی ناشی از بار ضربه انفجار و سپس تحلیل حرارتی ناشی از آتش سوزی

20- تحلیل اثرات خزش و جمع شدگی در قابهای بتنی

21- اندرکنش خاک و سازه با در نظر گرفتن پی

22- تحلیل استاتیکی و دینامیکی غیرخطی سازه بتنی با بتن فوق عملکرد بالا

23- تحلیل خرابی پیشرونده قاب بتنی

24- استفاده از روش‌های پردازش موازی برای اجرای هم‌زمان تحلیل‌های استاتیکی و دینامیکی غیرخطی سازه‌های بتنی

25- استفاده از روش‌های پردازش موازی برای اجرای هم‌زمان تحلیل بهینه‌یابی قطر میلگرد مقطع ستون بتنی (روش نیوتن- رافسون)

26- بهینه‌سازی ضریب محصورشدگی بتن در مقاطع مختلف ستون‌ها و تیرها با تحلیل استاتیکی غیرخطی و با هدف بیشینه‌سازی نسبت شکل‌پذیری (روش نیوتن- رافسون) و همینطور تحلیل حساسیت سازه با استفاده از ارتفاع مقطع ستون و ضریب محصور شدگی

27- تحلیل استاتیکی و دینامیکی غیرخطی سازه بتنی و میراگر ویسکوز به عنوان بادبند مهاری با استفاده از اوپنسیس و پایتون

28- مطالعه تطبیقی رفتار الاستیک و غیرالاستیک سازه‌ها با استفاده از تحلیل دینامیکی یا پوش‌آور با استفاده از اوپنسیس و پایتون

29- مطالعه تطبیقی رفتار الاستیک و غیرالاستیک سازه‌ها با استفاده از تحلیل دینامیکی تحت اثر بار ضربه و پوش‌آور با استفاده از اوپنسیس و پایتون

30- پیاده‌سازی تحلیل طیفی با ترکیب مودهای ارتعاشی به روش هم‌بستگی کامل برای قاب‌های بتن‌آرمه دوبعدی با استفاده از اوپنسیس

# CONCRETE FRAME EXAMPLES — OPENSEES_SALAR

This directory contains a collection of reinforced-concrete structural analysis examples implemented using **OpenSees/OpenSeesPy and Python** by **Salar Delavar Ghashghaei (Qashqai)**.

The following table maps the **30 analysis topics** documented in this project to their corresponding directories in the repository.

> **Note:** These examples are provided for educational and research purposes. The scripts may contain errors or inaccuracies and should be independently verified before use in engineering practice.

## Analysis Topics and Repository Links

| No. | Analysis Topic | Repository Link |
|---:|---|---|
| 1 | **Pushover Analysis** | [PUSHOVER](./PUSHOVER) |
| 2 | **Hysteretic Pushover Analysis** | [CYCLIC_PUSHOVER](./CYCLIC_PUSHOVER) |
| 3 | **Structural Analysis Considering Sectional Strengthening with Steel Plates or FRP Composites** | [COMPOSITE_SECTION (FRP or PLATE)](./COMPOSITE_SECTION%20(FRP%20or%20PLATE)) |
| 4 | **Thermal Analysis Due to Fire Exposure** | [THERMAL_LOAD](./THERMAL_LOAD) |
| 5 | **Seismic Response Analysis Under Earthquake Loading** | [SEISMIC](./SEISMIC) |
| 6 | **Incremental Dynamic Analysis (IDA)** | [SEISMIC_IDA](./SEISMIC_IDA) |
| 7 | **Computation of Acceleration, Velocity, and Displacement Response Spectra** | [SEISMIC_RESPONSE_SPECTRUM](./SEISMIC_RESPONSE_SPECTRUM) |
| 8 | **Free-Vibration / Modal Analysis** | [FREE-VIBRATION](./FREE-VIBRATION) |
| 9 | **Structural Analysis Under Blast, Harmonic, and Wind-Induced Impact Loads** | [IMPACT_LOAD](./IMPACT_LOAD) |
| 10 | **Post-Buckling Behavior Analysis of Frame Columns** | [POST_BUCKLING](./POST_BUCKLING) |
| 11 | **Sensitivity Analysis of Reinforcement Ratio and Section Dimensions on Structural Ductility Ratio and Structural Behavior Coefficient** | [SENSITIVITY](./SENSITIVITY) |
| 12 | **Optimization of Rebar Diameter and Section Depth Based on Structural Ductility Ratio, Structural Behavior Coefficient, and Structural Period** | [OPTIMIZATION](./OPTIMIZATION) |
| 13 | **Elastic Structural Analysis** | [ELASTIC_CONCRETE_FRAME](./ELASTIC_CONCRETE_FRAME) |
| 14 | **Structural Analysis Under Incrementally Increasing Distributed Loads** | [INCREMENTAL_DISTRIBUTED_LOAD](./INCREMENTAL_DISTRIBUTED_LOAD) |
| 15 | **Sequential Fire-Then-Earthquake Dynamic Analysis** | [THERMAL_&_SEISMSIC](./THERMAL_%26_SEISMSIC) |
| 16 | **Assessment of the Structural Ductility Damage Index** | [DUCTILITY_DAMAGE_INDEX](./DUCTILITY_DAMAGE_INDEX) |
| 17 | **Uncertainty Analysis of RC Frames and Probabilistic Seismic Assessment** | [UNCERTAINTY](./UNCERTAINTY) |
| 18 | **Evaluation of Structural Ductility Damage Index Considering Pinned Beam Connections** | [SIMPLY_SUPPORTED_BEAM](./SIMPLY_SUPPORTED_BEAM) |
| 19 | **Sequential Explosion Impact Load-Then-Fire Analysis** | [EXPLOSION_&_THERMAL](./EXPLOSION_%26_THERMAL) |
| 20 | **Analysis of Creep and Shrinkage Effects in Concrete Frames** | [CREEP_AND_SHRINKAGE](./CREEP_AND_SHRINKAGE) |
| 21 | **Soil–Structure Interaction with Foundation Consideration** | [SOIL-FOUNDATION](./SOIL-FOUNDATION) |
| 22 | **Nonlinear Static and Dynamic Analysis of an Ultra-High-Strength Concrete Frame** | [ULTRA_HIGH_STRENGTH_CONCRETE](./ULTRA_HIGH_STRENGTH_CONCRETE) |
| 23 | **Progressive Collapse of a Concrete Frame Structure** | [PROGRESSIVE_COLLAPSE](./PROGRESSIVE_COLLAPSE) |
| 24 | **Parallel Processing for Simultaneous Nonlinear Static and Dynamic Structural Analyses** | [PARALLEL_COMPUTING](./PARALLEL_COMPUTING) |
| 25 | **Parallel Processing for Simultaneous Optimization Analyses** | [PARALLEL_COMPUTING_&_OPTIMIZATION](./PARALLEL_COMPUTING_%26_OPTIMIZATION) |
| 26 | **Optimization of Concrete Confinement Enhancement Ratio Using Nonlinear Static Analysis and Newton–Raphson, with Sensitivity Analysis** | [CONFINEMENT_ENHANCEMENT_RATIO](./CONFINEMENT_ENHANCEMENT_RATIO) |
| 27 | **Nonlinear Static and Dynamic Analysis of an RC Frame with a Viscous Damper** | [VISCOUS_DAMPER](./VISCOUS_DAMPER) |
| 28 | **Comparative Study of Elastic and Inelastic Structural Behavior Using Dynamic or Pushover Analysis** | [ELASTIC_OR_INELASTIC_CONCRETE_FRAME](./ELASTIC_OR_INELASTIC_CONCRETE_FRAME) |
| 29 | **Comparative Study of Elastic and Inelastic Structural Behavior Under Impact Loading** | [ELASTIC_OR_INELASTIC_CONCRETE_FRAME_IMPACT_LOAD](./ELASTIC_OR_INELASTIC_CONCRETE_FRAME_IMPACT_LOAD) |
| 30 | **Response Spectrum Analysis with CQC Modal Combination** | [CQC_RESPONSE_SPECTRUM_ANALYSIS](./CQC_RESPONSE_SPECTRUM_ANALYSIS) |

## Supporting Directory

The following directory is a supporting resource rather than one of the 30 main analysis topics:

- [USEFUL_FUNCTIONS](./USEFUL_FUNCTIONS) — reusable Python/OpenSees functions used by different examples.

## Main Analysis Categories

The 30 examples cover the following major areas:

- **Nonlinear static analysis:** Pushover, cyclic pushover, incremental distributed loading, progressive collapse.
- **Dynamic and seismic analysis:** Earthquake response, IDA, response spectra, CQC modal combination, free vibration.
- **Fire and thermal analysis:** Fire exposure, sequential fire–earthquake analysis, explosion–fire analysis.
- **Extreme loading:** Blast/explosion, harmonic loading, wind-induced impact, and impact loading.
- **Material and section behavior:** FRP/steel-plate strengthening, confinement enhancement, ultra-high-strength concrete.
- **Structural performance:** Ductility, ductility damage index, elastic/inelastic behavior, post-buckling.
- **Uncertainty and reliability:** Probabilistic seismic assessment and uncertainty analysis.
- **Optimization and sensitivity:** Reinforcement/section optimization, confinement optimization, and parameter sensitivity.
- **Soil–structure interaction:** Foundation and soil–structure modeling.
- **Structural control:** Viscous damping systems.
- **High-performance computing:** Parallel processing for nonlinear analysis and optimization.

## Software and Methods

- Python
- OpenSees / OpenSeesPy
- Nonlinear Static Analysis
- Nonlinear Dynamic Analysis
- Pushover Analysis
- Cyclic Pushover Analysis
- Incremental Dynamic Analysis (IDA)
- Response Spectrum Analysis
- CQC Modal Combination
- Sensitivity Analysis
- Optimization
- Uncertainty Analysis
- Soil–Structure Interaction
- Structural Fire Analysis
- Impact and Blast Analysis
- Parallel Computing

## Author

**Salar Delavar Ghashghaei (Qashqai)**

Repository: [OPENSEES_SALAR](https://github.com/salardelavar/OPENSEES_SALAR)

This collection is intended as an educational and research-oriented resource for structural, civil, and earthquake engineering applications using OpenSees and Python.

