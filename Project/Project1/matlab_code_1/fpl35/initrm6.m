% Initialize rm6c model for Saab J35J Draken
global rm6 mrm6 arm6
mrm6(1,1:22)=0:0.1:2.1;
arm6(1,1:17)=0:1:16;
%load rm6.dat
rm6=[    39226.61544    39618.72596    40401.43414    41339.50941    42550.98121    44121.66749    45655.28018    46818.76780    48014.13137    49627.93101    51249.91070    52367.69128    52857.21703    52782.58627    52572.18998    52607.66748    52719.35899    52685.69026    52628.72842    52679.84787    52664.01080    52329.57632
    35026.08787    35761.82181    36710.05999    37699.11337    38883.14637    40392.85696    41939.22388    43253.12470    44564.75147    46115.03505    47688.54707    48967.63685    49682.12381    49716.34397    49525.22972    49557.20009    49664.02427    49631.03801    49574.49315    49620.44206    49603.22300    49284.24570
    30986.11536    31991.11970    33080.98662    34123.92119    35286.73443    36729.24153    38284.66848    39760.11729    41184.78330    42636.86671    44168.49389    45705.56977    46695.18731    46731.43118    46475.16485    46608.62182    46832.73680    46762.84472    46653.17993    46763.39399    46735.90801    46058.83544
    27782.18537    28845.70521    29903.74205    30948.81757    32083.88884    33404.88264    34869.16933    36401.21484    37932.42385    39427.92475    40980.80191    42594.92194    43789.32458    44170.05363    44175.99598    44306.14342    44467.06396    44489.01566    44489.44046    44587.99128    44625.98113    44375.13803
    24969.45410    25990.61766    26988.83166    28011.99287    29115.99747    30349.86879    31727.13952    33240.35955    34833.63823    36449.66143    38073.86154    39657.52905    40972.64027    41809.99449    42214.99706    42321.26162    42360.62889    42517.27197    42686.46637    42758.49722    42892.58564    43315.18674
    22014.92960    23027.32406    24090.63202    25103.89993    26178.11007    27435.66011    28832.67357    30306.96756    31889.40791    33598.91791    35313.60206    36893.13602    38264.29897    39353.16013    40015.84544    40192.33440    40227.96403    40449.49364    40704.01413    40818.50235    41019.14742    41631.94154
    19153.35843    20177.04836    21318.02727    22320.73592    23360.67941    24643.17620    26061.73402    27489.87773    29033.01164    30787.13763    32578.76842    34204.72087    35612.51746    36762.37209    37514.55878    37797.87890    37915.18321    38164.38017    38449.55968    38643.93044    38890.99480    39401.82854
    16697.58158    17709.81328    18808.48633    19780.61484    20776.70775    21975.89813    23308.32072    24687.53402    26209.78985    27958.39069    29782.14836    31487.77659    32948.08747    34064.63271    34787.82417    35135.41962    35345.70048    35641.44076    35962.86104    36224.42430    36520.11637    36988.80383
    14505.07910    15471.07971    16459.30981    17383.57378    18333.99248    19418.09346    20636.71463    21977.37809    23481.03245    25172.65172    26959.88489    28704.62276    30203.04822    31282.75226    31962.66609    32338.35114    32620.54930    32986.60410    33365.08507    33672.16944    34023.23919    34583.47747
    12370.79333    13256.24042    14144.02090    15008.44821    15917.75680    16949.38026    18123.62756    19444.68148    20909.34454    22504.60994    24185.61367    25869.68124    27354.75346    28454.28992    29163.21008    29560.96868    29883.70180    30328.42862    30779.01940    31116.84157    31526.40020    32267.98476
    10294.87667    11084.97749    11881.20135    12667.34385    13528.13210    14549.99208    15725.21376    17027.03767    18446.64250    19971.90870    21559.58569    23138.95621    24560.56732    25679.99273    26451.64849    26905.29971    27273.26534    27759.57235    28252.52601    28635.53421    29088.13720    29863.90830
     8348.54156     9044.73518     9744.50811    10427.67444    11213.76633    12215.38860    13377.71847    14629.57858    16001.80079    17519.58499    19083.59383    20570.81509    21888.07430    22961.87451    23767.59109    24322.43041    24762.05041    25217.87043    25685.90388    26139.99505    26600.71757    27100.32744
     6682.03049     7269.71691     7882.84026     8444.81520     9102.46543    10001.31853    11058.30857    12175.26782    13422.21502    14859.29661    16338.98336    17691.49589    18865.72983    19843.82844    20622.25035    21211.94942    21667.54478    22057.76687    22464.12577    22933.17051    23358.82607    23596.86144
     5372.05599     5824.73525     6366.16358     6799.81115     7270.37756     7948.43781     8766.83907     9637.58992    10627.07233    11793.59410    13008.79281    14126.60910    15116.86573    15969.98238    16642.88407    17111.64322    17462.41637    17792.84261    18136.40698    18500.02330    18848.47455    19136.01104
     4122.18059     4450.28002     4873.46861     5188.42491     5507.86178     5976.74224     6552.99592     7174.28119     7884.20491     8721.97834     9603.28157    10431.38463    11183.44840    11843.52350    12349.32856    12663.58739    12895.37651    13156.61747    13424.25930    13662.38249    13918.56591    14261.26304
     2681.71505     2920.75426     3150.22591     3365.21744     3628.47604     3995.84952     4427.92853     4880.35567     5384.23839     5966.96882     6565.61298     7108.07231     7575.91595     7965.27797     8276.88435     8516.69531     8707.01624     8873.86388     9041.75523     9224.93068     9398.02381     9525.76651
     1356.11302     1470.53652     1598.80140     1705.22551     1824.36316     1994.52481     2198.80910     2415.68655     2660.32704     2946.20898     3243.34557     3517.38225     3759.95854     3967.56851     4130.13027     4242.94250     4329.40208     4414.42603     4500.91376     4586.89818     4672.96396     4760.33371
    51896.11634    54869.27435    57871.85986    60787.69772    63716.50733    66753.58725    69760.65834    72732.51134    76429.79449    81542.72351    87713.92662    94120.08407    99121.67186   101460.49823   102219.90479   102706.11182   102775.32099   102117.68374   101209.16950   100481.80294    99406.01186    97211.82481
    45787.29494    48592.38767    51804.91923    55207.67413    58587.33759    61770.81746    64742.01152    67661.70080    71241.20969    76067.29689    81677.91513    87261.70099    91658.83250    94008.57925    94995.03575    95480.18673    95492.75299    94947.31861    94135.18487    93325.71778    92323.82499    90818.29955
    39732.79298    42971.86538    46458.01519    49825.80368    53244.80058    56882.88291    60366.14860    63449.73337    66946.70337    71569.44062    76569.66035    80965.85881    84304.32321    86403.75633    87644.73244    88395.75168    88415.14633    87536.17922    86502.00333    85885.76053    84670.65818    81442.42000
    34661.10853    37991.26069    41414.62955    44695.71987    48055.77149    51699.47001    55308.54906    58641.36769    62285.73875    66751.59635    71407.90603    75471.11294    78688.54936    80952.20975    82201.85029    82475.43363    82155.98678    81626.79597    80927.12006    80042.09072    79084.35816    78198.45224
    30270.89976    33430.87885    36654.34968    39850.04862    43105.92561    46508.21748    49957.09507    53402.48235    57179.38166    61551.24728    66110.26271    70334.14474    73914.01416    76567.69995    77906.45984    77799.14521    77251.55333    77155.34737    76808.70825    75596.48647    74843.69312    76360.37942
    26072.49151    29072.35215    32213.61946    35281.91734    38381.74259    41636.78257    44928.61364    48194.50630    51836.61850    56179.74845    60774.36805    65057.03795    68783.00013    71748.10072    73591.92101    74114.23838    73911.87951    73586.42149    72961.39362    71833.12471    70667.18944    70096.47396
    22052.71535    24905.86901    27960.20340    30847.46702    33754.11960    36894.94449    40043.30829    43039.95004    46436.51525    50708.90067    55319.10142    59605.35988    63424.80922    66720.15932    69257.53628    70780.53677    71119.22368    70304.59638    69084.93411    68044.15898    66391.48480    62992.44830
    18253.29815    20935.49336    23754.85679    26424.89241    29130.31612    32072.57203    35048.80504    37914.05441    41159.23599    45223.98439    49706.93289    54088.28054    58215.49093    61980.81968    65088.42582    67177.60706    67816.31222    66841.18468    65234.99483    63903.23274    62296.14040    59500.14773
    14459.04149    17012.23496    19598.97758    22117.58962    24678.62353    27397.77033    30199.04315    33041.90014    36219.25736    40003.75275    44253.45379    48725.89700    53191.06548    57378.93908    60847.03822    63151.51639    64015.51891    63394.44744    62005.73714    60601.97229    59314.30017    58118.50883
    10563.36517    13085.11645    15598.71023    18105.11945    20621.67681    23178.13761    25837.58774    28677.59640    31800.33552    35306.97762    39270.09473    43684.06064    48259.05668    52604.39878    56210.13394    58642.38958    59890.88368    60109.62082    59686.16143    58966.37517    57895.81277    56319.94526
     7151.08739     9648.26603    12138.54262    14646.50008    17138.45365    19601.00812    22160.19402    24941.52631    27929.03187    31119.03306    34698.52226    38800.71282    43153.03213    47359.26259    50934.39278    53500.08166    55193.46363    56252.89291    56806.12765    56841.94109    55903.76307    53423.68756
     4804.16084     7182.75605     9555.79733    11942.74690    14336.91931    16740.95942    19220.98080    21835.14832    24546.36344    27343.72126    30416.35512    33924.63770    37713.86095    41518.08527    44941.52553    47648.18751    49671.08628    51107.11864    51939.69921    52081.85691    51278.56002    49232.76164
     3401.67797     5572.32973     7720.79229     9883.40674    12084.64671    14341.80659    16655.33127    19017.94818    21408.36436    23831.95456    26414.78481    29275.90067    32385.57464    35633.20267    38733.44765    41395.53128    43481.64366    44900.86826    45596.89381    45550.55595    44856.67276    43638.55794
     2624.90709     4533.56255     6386.33771     8280.51514    10206.55444    12143.00107    14127.56678    16196.46131    18320.72000    20475.72784    22719.44272    25116.16633    27673.00248    30345.39592    32939.35301    35234.25661    37052.43398    38256.21604    38825.00375    38771.45854    38124.21357    36915.89494
     2066.68011     3711.68800     5302.21757     6947.53090     8597.33615    10197.30249    11836.49704    13607.39472    15472.70363    17384.25594    19380.14712    21502.21487    23721.00267    25971.44621    28117.34433    30004.55578    31478.31577    32423.53700    32884.46566    32892.68673    32269.80700    30784.93869
     1501.05347     2918.84331     4337.68685     5787.23741     7223.86350     8612.00273    10021.65376    11525.26572    13099.52835    14723.80823    16483.93854    18450.65893    20527.86814    22582.43737    24499.96854    26164.01182    27431.17968    28196.68185    28537.05479    28513.18359    27942.02107    26579.53697
     1515.41796     2706.12241     3842.91084     5058.29299     6271.94499     7399.21443     8550.96830     9845.66122    11236.57896    12671.65216    14242.55979    16022.32535    17865.60194    19590.66387    21088.64057    22291.72195    23223.48472    23915.10339    24336.75678    24410.47454    23926.68518    22642.91690
    52656.40890    57446.91201    61462.24338    64946.36187    68143.22636    71296.79573    74651.02884    78449.88456    82937.32177    88134.54841    93171.76873    96157.91223    92007.81296    74837.78115    38764.12702   -22096.83922  -113628.80733  -241715.46710  -412240.50832  -631087.62076  -904140.49419 -1237282.81841
    42853.40373    49140.82927    54269.69907    58549.08163    62288.04545    65795.65904    69380.99089    73353.10952    78021.08343    83559.67992    89606.46157    95016.27694    96050.32337    88321.38542    67442.24764    29025.69458   -31315.48922  -117968.51920  -235320.61083  -387758.97954  -579670.84079  -815443.41003
    35797.96091    42571.62827    48071.60672    52625.26017    56559.95250    60203.04762    63881.90941    67923.90178    72656.38862    78330.50177    84892.44481    91703.92255    96093.57298    94881.76734    84888.87689    62935.27289    25841.32660   -29572.59074  -106486.10786  -208078.85350  -337530.45642  -498020.54535
    30767.43221    37307.06382    42645.54464    47095.66834    50970.22861    54582.01911    58243.83351    62268.46550    66968.70873    72614.02795    79300.57220    86700.69238    92964.86264    95863.08796    93165.37331    82641.72364    62062.14394    29196.63917   -18184.78571   -82312.12572  -165415.37590  -269724.53127
    27039.16943    32914.89077    37769.09115    41881.07704    45530.15489    48995.63116    52556.81231    56493.00481    61083.51511    66577.27249    73101.69751    80486.42970    87491.49320    92609.50834    94333.09542    91154.87475    81567.46663    64063.49139    37135.56933     -723.67925   -51021.63403  -115265.67469
    23890.52433    28962.86396    33219.82459    36902.25712    40251.01243    43506.94141    46910.89492    50703.72387    55126.27914    60387.24939    66566.67450    73540.97780    80500.76552    86465.18954    90453.40177    91484.55412    88577.79850    80752.28681    67027.17095    46421.60284    17954.73438   -19354.28253
    20598.84869    25018.73821    28775.32328    32079.97947    35144.08235    38179.00749    41396.13047    45006.82684    49222.47218    54210.97267    59966.35691    66344.17999    72819.98043    78774.29260    83587.65089    86640.58968    87313.64335    84987.34629    79042.23287    68858.83748    53817.69450    33299.33831
    16441.49431    20650.26837    24213.16553    27335.01494    30220.64574    33074.88706    36102.56805    39508.51786    43497.56561    48215.45633    53571.59852    59375.87955    65276.43879    70880.97859    75797.20132    79632.80934    81995.50502    82492.99072    80732.96879    76323.14161    68871.21153    57984.88092
    11096.94973    15672.32268    19450.50591    22658.27556    25522.40777    28269.67868    31126.86443    34320.74114    38078.08496    42552.21887    47602.65323    53006.47317    58504.87549    63830.08494    68714.32625    72889.82414    76088.80337    78043.48865    78486.10473    77148.87634    73764.02821    68063.78508
     5848.25066    10888.22300    14962.80383    18321.23789    21212.76992    23886.64465    26592.10681    29578.40115    33094.77240    37310.79666    42077.37554    47168.57108    52369.76162    57469.15460    62254.95744    66515.37757    70038.62241    72612.89939    74026.41594    74067.37948    72523.99744    69184.47725
     2379.56958     7348.40461    11365.09495    14665.51964    17485.55770    20061.08816    22627.99006    25422.14241    28679.42425    32565.23055    36963.02009    41685.33690    46543.00229    51346.40689    55905.94132    60031.99622    63534.96222    66225.22993    67913.19000    68409.23304    67523.74970    65067.13058
     1818.10878     5763.47196     9079.59150    11929.94161    14477.99650    16887.23037    19321.11743    21943.13189    24916.74796    28363.25214    32235.18003    36429.09838    40786.44923    45134.89344    49302.09188    53115.70543    56403.39494    58992.82128    60711.64532    61387.52794    60848.12999    58921.11235
     3063.18980     5484.70625     7757.21201     9942.13698    12100.91102    14294.96402    16585.72584    19034.62636    21703.09545    24647.35595    27902.80249    31468.83973    35221.74062    39006.99517    42670.09342    46056.52538    49011.78107    51381.35052    53010.72374    53745.39077    53430.84162    51912.56630
     4457.16399     5523.55782     6856.05159     8426.94201    10208.52584    12173.09979    14292.96060    16540.40501    18887.72974    21333.72728    23983.17310    26922.70926    30060.46123    33259.92530    36384.59778    39297.97499    41863.55321    43944.82876    45405.29794    46108.45705    45917.80241    44696.83031
     4342.38268     4891.47705     5834.20533     7109.19298     8655.06549    10410.44835    12313.96704    14304.24705    16319.91385    18338.55139    20493.57747    22908.85529    25514.19584    28190.89705    30820.25685    33283.57317    35462.14394    37237.26710    38490.24058    39102.36231    38954.93021    37929.24224
     1061.19722     2599.91431     4149.76834     5713.72615     7294.75455     8895.82039    10519.89049    12169.93169    13848.91081    15578.01357    17451.30124    19545.42615    21794.52924    24097.12364    26351.72247    28456.83883    30310.98586    31812.67666    32860.42435    33352.74205    33188.14288    32265.13996
    -7044.04103    -2339.68004     1260.83574     3965.37778     5981.81757     7518.02657     8781.87627     9981.23813    11323.98364    12968.29909    14873.63004    16950.57015    19113.04621    21275.81831    23353.64651    25261.29090    26913.51155    28225.06853    29110.72192    29485.23181    29263.35826    28359.86135
        1.26114        1.32216        1.38911        1.45552        1.52468        1.60010        1.67649        1.74985        1.83032        1.92800        2.03887        2.15566        2.27226        2.38469        2.49626        2.61065        2.72584        2.83933        2.95257        3.06708        3.18079        3.29076
        1.11978        1.18052        1.24300        1.30614        1.37197        1.44229        1.51498        1.58852        1.66757        1.75763        1.86144        1.97787        2.09308        2.19502        2.29141        2.39189        2.49384        2.59318        2.69203        2.79277        2.89236        2.98643
        0.98355        1.04302        1.10044        1.16051        1.22350        1.28894        1.35799        1.43172        1.50927        1.59152        1.68822        1.80452        1.91821        2.00883        2.09040        2.17980        2.27101        2.35576        2.43972        2.52860        2.61437        2.68556
        0.86216        0.92020        0.97426        1.03112        1.09111        1.15366        1.22089        1.29440        1.37068        1.44788        1.53561        1.64016        1.74315        1.82664        1.89944        1.97375        2.04880        2.12166        2.19391        2.26734        2.34012        2.40952
        0.75674        0.81090        0.86197        0.91598        0.97367        1.03485        1.10106        1.17319        1.24781        1.32233        1.40195        1.49064        1.57989        1.66014        1.73016        1.79195        1.85211        1.91601        1.97975        2.03964        2.10213        2.17620
        0.66646        0.71261        0.75915        0.81077        0.86720        0.92737        0.99189        1.06117        1.13325        1.20630        1.28145        1.35988        1.44005        1.51887        1.58973        1.64822        1.70205        1.75921        1.81659        1.87057        1.92675        1.99302
        0.59089        0.62884        0.66898        0.71618        0.76950        0.82707        0.88896        0.95541        1.02525        1.09708        1.16983        1.24290        1.31726        1.39261        1.46206        1.51984        1.57113        1.62235        1.67394        1.72512        1.77647        1.82889
        0.52892        0.56232        0.59490        0.63384        0.67966        0.73152        0.78968        0.85415        0.92268        0.99279        1.06316        1.13295        1.20189        1.26910        1.33052        1.38293        1.42957        1.47467        1.51999        1.56616        1.61171        1.65464
        0.47605        0.50647        0.53387        0.56405        0.60017        0.64429        0.69659        0.75642        0.82143        0.88907        0.95755        1.02508        1.08920        1.14768        1.19970        1.24536        1.28705        1.32735        1.36752        1.40816        1.44850        1.48743
        0.42667        0.45320        0.48041        0.50513        0.53273        0.56902        0.61312        0.66315        0.72000        0.78411        0.85134        0.91684        0.97744        1.03107        1.07828        1.12046        1.15966        1.19778        1.23535        1.27268        1.31022        1.34849
        0.37405        0.39774        0.42577        0.44883        0.47242        0.50360        0.54091        0.58213        0.63049        0.68872        0.75190        0.81396        0.87184        0.92376        0.97000        1.01150        1.04970        1.08610        1.12193        1.15806        1.19405        1.22911
        0.31385        0.33754        0.36313        0.38725        0.41281        0.44329        0.47837        0.51725        0.56091        0.61028        0.66397        0.71981        0.77466        0.82564        0.87199        0.91356        0.95037        0.98330        1.01627        1.05210        1.08608        1.11157
        0.25241        0.27579        0.29688        0.32070        0.34802        0.37863        0.41267        0.45006        0.48940        0.52970        0.57280        0.62020        0.66915        0.71626        0.75986        0.79881        0.83232        0.86092        0.89006        0.92361        0.95426        0.97190
        0.19683        0.21644        0.23314        0.25216        0.27526        0.30283        0.33292        0.36362        0.39527        0.42868        0.46428        0.50206        0.54071        0.57859        0.61394        0.64530        0.67276        0.69729        0.72205        0.74923        0.77483        0.79332
        0.14671        0.16088        0.17377        0.18721        0.20363        0.22458        0.24737        0.26937        0.29234        0.31836        0.34624        0.37432        0.40218        0.42949        0.45493        0.47741        0.49779        0.51726        0.53654        0.55601        0.57538        0.59415
        0.09922        0.10859        0.11795        0.12839        0.13982        0.15201        0.16549        0.18066        0.19660        0.21255        0.22954        0.24853        0.26834        0.28744        0.30507        0.32080        0.33471        0.34723        0.35965        0.37289        0.38562        0.39592
        0.04926        0.05398        0.05845        0.06329        0.06890        0.07546        0.08265        0.09013        0.09795        0.10627        0.11516        0.12459        0.13420        0.14354        0.15220        0.15990        0.16677        0.17313        0.17943        0.18601        0.19242        0.19803
        3.26673        3.45804        3.67325        3.88591        4.09876        4.32129        4.56095        4.82528        5.12425        5.46396        5.83247        6.21449        6.59912        6.97912        7.35742        7.73806        8.11948        8.49959        8.87973        9.26120        9.64167       10.01790
        2.90238        3.13347        3.37505        3.60778        3.83061        4.04740        4.26346        4.48881        4.75079        5.07295        5.44616        5.84811        6.23647        6.57817        6.89734        7.22446        7.55832        7.89121        8.22113        8.54868        8.87930        9.21964
        2.66644        2.86996        3.08717        3.29398        3.50810        3.74746        3.98893        4.21205        4.44830        4.73204        5.05726        5.40702        5.76087        6.10234        6.43419        6.76154        7.07950        7.38549        7.69609        8.02249        8.33511        8.59419
        2.44158        2.61947        2.80235        2.98187        3.17508        3.39755        3.63351        3.86675        4.11068        4.38076        4.67122        4.97384        5.29204        5.62483        5.94203        6.22019        6.49187        6.78959        7.08967        7.36907        7.66362        8.02384
        2.20204        2.36055        2.51270        2.66979        2.83867        3.02487        3.23299        3.46603        3.72150        3.99433        4.27453        4.55554        4.84937        5.16042        5.44391        5.66684        5.89237        6.17774        6.46041        6.68590        6.95204        7.39467
        1.94891        2.08540        2.21539        2.35195        2.50011        2.66327        2.84619        3.05428        3.29398        3.56707        3.85560        4.14150        4.42600        4.70934        4.96863        5.18850        5.40676        5.65749        5.90649        6.12197        6.35004        6.65626
        1.69927        1.81437        1.92851        2.04527        2.17386        2.32182        2.48437        2.65970        2.87020        3.13414        3.42925        3.72647        4.01019        4.27205        4.51902        4.75899        4.98815        5.20276        5.41131        5.62134        5.82452        6.00855
        1.46687        1.56652        1.66815        1.76636        1.87401        2.00359        2.14798        2.30230        2.49037        2.73365        3.01576        3.31109        3.59523        3.85166        4.09423        4.33722        4.56628        4.76523        4.93940        5.10668        5.31395        5.61531
        1.23647        1.32750        1.42139        1.50791        1.60007        1.71232        1.84178        1.98559        2.15838        2.37460        2.63053        2.91424        3.19919        3.46179        3.70493        3.93555        4.15019        4.34122        4.49472        4.61999        4.82570        5.24531
        0.99884        1.08700        1.17891        1.26432        1.35219        1.45410        1.57298        1.71030        1.87035        2.05906        2.28595        2.55424        2.83815        3.10878        3.35368        3.56776        3.75883        3.93380        4.08307        4.20860        4.37517        4.66325
        0.79804        0.88412        0.97472        1.06076        1.14768        1.24383        1.35438        1.48266        1.62756        1.79048        1.98712        2.22884        2.49531        2.75948        2.99925        3.19992        3.37150        3.52840        3.67812        3.82163        3.94073        4.01245
        0.67667        0.75687        0.84157        0.92105        1.00167        1.09202        1.19339        1.30558        1.42969        1.56892        1.73372        1.93202        2.15456        2.38646        2.60759        2.80084        2.96656        3.11197        3.25438        3.39946        3.49643        3.48037
        0.61835        0.68921        0.76377        0.83291        0.90323        0.98284        1.07032        1.16327        1.26485        1.37953        1.51126        1.66316        1.83538        2.02467        2.21718        2.39794        2.55831        2.69524        2.82181        2.94557        3.03586        3.05246
        0.59214        0.65192        0.71399        0.77385        0.83316        0.89507        0.96248        1.03746        1.11900        1.20682        1.30659        1.42439        1.56188        1.71718        1.87842        2.03296        2.17471        2.29958        2.40469        2.49090        2.57283        2.66849
        0.56804        0.61691        0.66677        0.71669        0.76529        0.81228        0.86203        0.91870        0.98092        1.04780        1.12614        1.22326        1.34108        1.47715        1.61702        1.74628        1.86257        1.96532        2.04899        2.11310        2.18233        2.28766
        0.53099        0.56969        0.60855        0.64582        0.68377        0.72444        0.76496        0.80363        0.84827        0.90710        0.98028        1.06749        1.17446        1.30330        1.43561        1.55135        1.64439        1.71600        1.78335        1.85791        1.91272        1.91116
        0.52491        0.54889        0.57139        0.59654        0.62110        0.64244        0.66769        0.70429        0.75140        0.80734        0.87565        0.95988        1.05847        1.16636        1.26948        1.35532        1.42661        1.48905        1.54525        1.59680        1.64448        1.68883
        3.32283        3.62416        3.90285        4.16876        4.43175        4.70171        4.98848        5.30193        5.65193        6.02702        6.33049        6.46363        6.40502        6.15260        5.70427        5.05798        4.21163        3.16315        1.91046        0.45149       -1.21584       -3.09362
        3.09993        3.32575        3.56149        3.80957        4.07244        4.35250        4.65219        4.97393        5.32016        5.68390        6.02067        6.28055        6.42958        6.43776        6.27509        5.91159        5.31725        4.46210        3.31614        1.84937        0.03181       -2.16654
        2.83686        3.01497        3.21752        3.44374        3.69288        3.96418        4.25689        4.57025        4.90351        5.25375        5.60942        5.95154        6.24009        6.42977        6.47530        6.33139        5.95275        5.29409        4.31013        2.95557        1.18513       -1.04648
        2.54655        2.69839        2.87570        3.07767        3.30353        3.55246        3.82369        4.11642        4.42985        4.76457        5.12664        5.51385        5.88528        6.19036        6.37850        6.39912        6.20164        5.73548        4.95005        3.79479        2.21909        0.17239
        2.24194        2.38256        2.54078        2.71780        2.91483        3.13306        3.37370        3.63795        3.92701        4.24432        4.60220        5.00472        5.41389        5.78123        6.05826        6.19649        6.14744        5.86261        5.29353        4.39171        3.10866        1.39590
        1.93595        2.07403        2.21752        2.37054        2.53722        2.72167        2.92802        3.16037        3.42285        3.72099        4.06598        4.46139        4.87463        5.26410        5.58817        5.80523        5.87366        5.75185        5.39818        4.77104        3.82881        2.52987
        1.64152        1.77935        1.91067        2.04232        2.18116        2.33402        2.50775        2.70920        2.94521        3.22257        3.54787        3.92112        4.31625        4.70067        5.04180        5.30705        5.46383        5.47955        5.32163        4.95748        4.35451        3.48013
        1.37157        1.50508        1.62499        1.73957        1.85708        1.98581        2.13401        2.30997        2.52194        2.77704        3.07775        3.42116        3.78747        4.15267        4.49275        4.78369        5.00147        5.12208        5.12150        4.97572        4.66072        4.15249
        1.13590        1.25733        1.36576        1.46889        1.57439        1.68993        1.82319        1.98186        2.17360        2.40493        2.67774        2.98963        3.32576        3.66817        3.99892        4.30004        4.55360        4.74164        4.84619        4.84932        4.73307        4.47949
        0.93180        1.04044        1.14044        1.23774        1.33828        1.44800        1.57283        1.71870        1.89156        2.09696        2.33895        2.61822        2.92355        3.24073        3.55556        3.85384        4.12136        4.34391        4.50730        4.59733        4.59977        4.50045
        0.75343        0.85831        0.95698        1.05375        1.15292        1.25881        1.37574        1.50800        1.65991        1.83639        2.04472        2.28950        2.56203        2.85029        3.14228        3.42598        3.68937        3.92045        4.10721        4.23763        4.29969        4.28140
        0.59795        0.71442        0.82100        0.92146        1.01958        1.11914        1.22392        1.33771        1.46429        1.60872        1.78117        1.98941        2.22642        2.48152        2.74400        3.00316        3.24831        3.46875        3.65378        3.79270        3.87481        3.88943
        0.47448        0.61050        0.72862        0.83304        0.92796        1.01757        1.10607        1.19765        1.29651        1.40841        1.54545        1.71737        1.91820        2.13806        2.36706        2.59532        2.81296        3.01009        3.17683        3.30328        3.37958        3.39583
        0.39518        0.54783        0.67359        0.77758        0.86491        0.94070        1.01006        1.07810        1.14993        1.23217        1.33751        1.47620        1.64289        1.82830        2.02316        2.21819        2.40414        2.57172        2.71166        2.81469        2.87155        2.87295
        0.37218        0.52771        0.64966        0.74415        0.81726        0.87509        0.92375        0.96932        1.01791        1.07670        1.15726        1.26871        1.40601        1.56062        1.72397        1.88751        2.04266        2.18087        2.29359        2.37224        2.40828        2.39314
        0.41761        0.55143        0.65060        0.72183        0.77184        0.80732        0.83501        0.86160        0.89381        0.93869        1.00464        1.09769        1.21311        1.34343        1.48120        1.61899        1.74934        1.86479        1.95791        2.02124        2.04734        2.02875
        0.54361        0.62030        0.67015        0.69970        0.71547        0.72396        0.73170        0.74521        0.77099        0.81484        0.87957        0.96597        1.06969        1.18509        1.30653        1.42838        1.54499        1.65072        1.73994        1.80701        1.84628        1.85213];