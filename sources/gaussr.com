 %NprocShared=2                                                                                                    
 %Mem=1000MB                                                                                                      
 %Chk=AzDCM.chk                                                                                                    
#p blyp sddall scf=tight force NoSymm units=au                                                                     
                                                                                                                   
AZ-S16 anion NO_BROKEN input optimised pm6 geometry empiricaldispersion=gd3                                        
                                                                                                                   
 -1  2                                                                                                             
    7  -5.871312957679	-0.098326230200	4.810437713630
    6  4.801119474064	-5.045532872773	6.808414919407
    6  2.373554605715	-4.466826920973	7.603286643666
    6  0.435508174550	-3.219955382194	6.320487631629
    6  5.960164101043	-4.568225845867	4.416388240733
    6  0.495036437494	-2.236883942532	3.831058798834
    6  4.987612766769	-3.390359546594	2.350137119464
    6  2.539289256839	-2.266312647615	2.068739890588
    6  -1.596937635910	-0.912637567199	2.691126861973
    6  -3.916529222145	-0.462092841802	3.857458272925
    6  -0.825999293137	-0.062287263100	0.230045810913
    6  1.731782823648	-0.904872682515	-0.143995241679
    6  3.045579013460	-0.319677411171	-2.400524777099
    6  1.852067671524	1.055211734825	-4.262422803972
    6  -0.699965898361	1.916152064178	-3.957292504580
    6  -2.075499440969	1.366801118143	-1.686597582060
    6  -4.603146095096	2.269704706026	-1.465938040854
    6  -1.909696759705	3.313500605275	-5.918174386348
    6  -5.714078294706	3.662733333181	-3.380761627560
    6  -4.366877416024	4.162851244249	-5.640252363832
    1  6.011250957348	-6.025216033477	8.168713489923
    1  1.893545270438	-5.067138111662	9.537878655459
    1  -1.316126222135	-2.907460820058	7.378824973479
    1  7.905848805242	-5.288623463189	4.274636103975
    1  6.188195463898	-3.250948780516	0.662144919517
    1  4.994062402064	-0.938870745390	-2.628272681033
    1  2.842138656787	1.529527325261	-6.001063443760
    1  -5.663511113088	1.850125033059	0.244555128169
    1  -0.870606278528	3.689559885378	-7.656017561707
    1  -7.620808184380	4.413213500828	-3.165684337354
    1  -5.294003513492	5.201542993211	-7.152632314153
   17  10.878170695395	-6.627252544127	0.067587944906
   17  9.740085020967	-11.904257971044	-0.116178472987
    6  8.439123851056	-8.856796790354	-0.531750036824
    1  7.773142787488	-8.634047212044	-2.507039190637
    1  6.867903498103	-8.560642690098	0.833870002485
   17  -11.224808829319	4.999965906506	1.416457451764
   17  -8.990443891723	8.179089976124	5.182596598993
    6  -8.421025943872	5.587827581425	3.140902468834
    1  -6.876418604063	6.075853134962	1.803454466729
    1  -7.936099541194	3.891707669376	4.289540373993
   17  4.527997355928	-6.026553959580	-10.538747535851
   17  6.972420471118	-4.653390133545	-5.936425361349
    6  3.996795339683	-4.758200014110	-7.493846934073
    1  2.686368531568	-5.994740420116	-6.416960040492
    1  3.223814202960	-2.812126026328	-7.642409643817
   17  2.297111404146	6.541523228335	-4.257142909154
   17  4.601131647037	9.406841265881	-8.216944970039
    6  4.074795675587	9.290035403817	-4.912736148164
    1  5.919565224782	9.191843344172	-3.920529772370
    1  3.000728253398	10.992809041583	-4.327247969270
   17  2.614494687790	-9.556103174295	-2.480841917242
   17  2.144215552376	-9.916504183652	2.886961070957
    6  0.871057923074	-8.451000779808	0.171469969940
    1  -1.134549996831	-9.009265564017	-0.057810501889
    1  1.076410682595	-6.356719350794	0.323708197015
   17  4.913034724893	12.309887685690	0.444405005188
   17  7.813038685627	7.813531904148	1.006093973152
    6  4.857219773802	9.212564191231	1.722300177908
    1  3.310463815378	8.098481818156	0.837776066404
    1  4.595694904950	9.308352519233	3.799882097076
   17  -2.014877026592	-3.786468820978	-8.406161358102
   17  -6.425794510168	-1.060204391271	-9.887610271064
    6  -4.887015085755	-2.477175848543	-7.270842242735
    1  -4.462484330318	-1.014760257202	-5.819617609559
    1  -6.091857225590	-3.999896380269	-6.482233070883
   17  -1.861450161779	12.758537015775	-3.353851927413
   17  -6.168959941580	9.826330231993	-1.924430954371
    6  -2.923017045593	9.645421080014	-2.700877848824
    1  -1.849490085078	8.860319239868	-1.078495050521
    1  -2.670666797121	8.446017571993	-4.406060887411
   17  -3.988586369357	-7.368376677685	9.150908096659
   17  -5.292986840832	-10.385400039535	4.864026567336
    6  -4.157476863124	-7.371963377887	5.813300255093
    1  -2.238905949495	-7.063994709842	5.005162653374
    1  -5.494373065198	-5.883977902228	5.174852501296
   17  -11.503237298438	-4.799282660262	-5.408632411049
   17  -9.121141130598	-0.590293752453	-2.999161670433
    6  -12.022450892090	-2.194136447656	-3.372142586663
    1  -12.702070109167	-2.890311555393	-1.512901514734
    1  -13.430370561181	-0.908642686152	-4.234980201055
   17  -1.723299843033	4.938317370839	2.513929132121
   17  -0.804722866595	10.054920953783	3.985489108224
    6  -1.873610549451	7.065279723613	5.070845754337
    1  -0.629484782744	6.413522628654	6.629093137381
    1  -3.846324006544	7.248526466819	5.750922285139
   17  0.615969461435	2.116669014522	9.339200408643
   17  4.508202474676	5.240628198839	7.303659226776
    6  2.952768893845	2.268837321730	6.955841903745
    1  2.055846289176	2.158231651112	5.051665034078
    1  4.337337217099	0.708981781746	7.202111013518
   17  8.600534808785	0.059184332788	4.142812588325
   17  13.088851695354	2.497745517510	5.940384337600
    6  10.357347055077	2.895315550195	4.042545609385
    1  9.201295754294	4.460605942541	4.821197814237
    1  10.934012102652	3.296084889225	2.066336158945
   17  -14.445862142397	-1.367481419551	3.322352082487
   17  -10.666459007433	-5.188929071267	2.996389442194
    6  -11.161788352479	-1.939185936024	3.674500657817
    1  -10.552408366073	-1.527664055293	5.641671548158
    1  -10.077265075590	-0.752233833499	2.323560011116
   17  -4.989668788802	-8.460216974163	-3.543034300406
   17  -4.000778882374	-4.631297345313	0.158721877440
    6  -6.306834360311	-6.745638217238	-0.979199390814
    1  -6.868415613885	-8.092463040419	0.521422683502
    1  -7.975024120101	-5.675637485690	-1.656471568033
   17  14.256991574216	0.316729438402	-1.389891681773
   17  9.506734079046	2.010593017447	-3.307642454260
    6  11.167385721287	-0.671604888543	-2.210597852009
    1  10.206764228096	-1.445273885580	-0.510815650713
    1  11.268051432441	-2.115758166526	-3.723935893480
       