#!/bin/csh
# ======================================================================= #
# Copyright 1996-2017 the Authors
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# ======================================================================= #

# script that writes a ferret script to plot the wanted reaction rates and executes it

# name of plot
# ------------
set pl_name = ("rxn_atm")

# define number of runs and directories of netCDF files
# -----------------------------------------------------

# max number of runs is 1 because of number of variables in rxnrate!
set n_runs = ("1")
set m_runs = ("1")
set runs = ("/data/Elise/Mistra_Elise_v008.4/")


# how to plot
# -----------

# 1 - evolution w/ time; 2 - vertical profile; 3 - contour plot; 4 - box model
 set pl_ty = ("1")
# set pl_ty = ("2")
#set pl_ty = ("3")
# set pl_ty = ("4")

# times (h from model start)/altitudes (m); n_pl: number of heights/times to be overplotted
set n_pl   = ("2")
set m_pl   = ("1" "2")

# so far put levels and timesteps instead of height and time...
# set pld = ("240" "288")
# set pld = ("2" "6" "8")
# set pld = ("192" "240")
set pld = ("1" "5")

# maximum and minimum heights if plot type = 2,3 [now in levels, later in m]
# set minheight = 2 - doesn't work!! set to 2 by hand
set maxheight = 15

# NOTE: if you chose more than 4 lines it'll be hard to tell the differences!

# unit of plot: 1 - mol m-3 s-1; 2 - mol mol-1 s-1; 3 - molec cm-3 s-1
# --------------------------------------------------------
set pl_un = 1


# what to plot
# ------------

set n_spec = ('1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' '23' '24' '25' '26' '27' '28' '29' '30' '31' '32' '33' '34' '35' '36' '37' '38' '39' '40' '41' '42' '43' '44' '45' '46' '47' '48' '49' '50' '51' '52' '53' '54' '55' '56' '57' '58' '59' '60' '61' '62' '63' '64' '65' '66' '67' '68' '69' '70' '71' '72' '73' '74' '75' '76' '77' '78' '79' '80' '81' '82' '83' '84' '85' '86' '87' '88' '89' '90' '91' '92' '93' '94' '95' '96' '97' '98' '99' '100' '101' '102' '103' '104' '105' '106' '107' '108' '109' '110' '111' '112' '113' '114' '115' '116' '117' '118' '119' '120' '121' '122' '123' '124' '125' '126' '127' '128' '129' '130' '131' '132' '133' '134' '135' '136' '137' '138' '139' '140' '141' '142' '143' '144' '145' '146' '147' '148' '149' '150' '151' '152' '153' '154' '155' '156' '157' '158' '159' '160' '161' '162' '163' '164' '165' '166' '167' '168' '169' '170' '171' '172' '173' '174' '175' '176' '177' '178' '179' '180' '181' '182' '183' '184' '185' '186' '187' '188' '189' '190' '191' '192' '193' '194' '195' '196' '197' '198' '199' '200' '201' '202' '203' '204' '205' '206' '207' '208' '209' '210' '211' '212' '213' '214' '215' '216' '217' '218' '219' '220' '221' '222' '223' '224' '225' '226' '227' '228' '229' '230' '231' '232' '233' '234' '235' '236' '237' '238' '239' '240' '241' '242' '243' '244' '245' '246' '247' '248' '249' '250' '251' '252' '253' '254' '255' '256' '257' '258' '259' '260' '261' '262' '263' '264' '265' '266' '267' '268' '269' '270' '271' '272' '273' '274' '275' '276' '277' '278' '279' '280' '281' '282' '283' '284' '285' '286' '287' '288' '289' '290' '291' '292' '293' '294' '295' '296' '297' '298' '299' '300' '301' '302' '303' '304' '305' '306' '307' '308' '309' '310' '311' '312' '313' '314' '315' '316' '317' '318' '319' '320' '321' '322' '323' '324' '325' '326' '327' '328' '329' '330' '331' '332' '333' '334' '335' '336' '337' '338' '339' '340' '341' '342' '343' '344' '345' '346' '347' '348' '349' '350' '351' '352' '353' '354' '355' '356' '357' '358' '359' '360' '361' '362' '363' '364' '365' '366' '367' '368' '369' '370' '371' '372' '373' '374' '375' '376' '377' '378' '379' '380' '381' '382' '383' '384' '385' '386' '387' '388' '389' '390' '391' '392' '393' '394' '395' '396' '397' '398' '399' '400' '401' '402' '403' '404' '405' '406' '407' '408' '409' '410' '411' '412' '413' '414' '415' '416' '417' '418' '419' '420' '421' '422' '423' '424' '425' '426' '427' '428' '429' '430' '431' '432' '433' '434' '435' '436' '437' '438' '439' '440')
#'441' '442' '443' '444' '445' '446' '447' '448' '449' '450' '451' '452' '453' '454' '455' '456' '457' '458' '459' '460' '461' '462' '463' '464' '465' '466' '467' '468' '469' '470' '471' '472' '473' '474' '475' '476' '477' '478' '479' '480' '481' '482' '483' '484' '485' '486' '487' '488' '489' '490' '491' '492' '493' '494' '495' '496' '497' '498' '499' '500' '501' '502' '503' '504' '505' '506' '507' '508' '509' '510' '511' '512' '513' '514' '515' '516' '517' '518' '519' '520' '521' '522' '523' '524' '525' '526' '527' '528' '529' '530' '531' '532' '533' '534' '535' '536' '537' '538' '539' '540' '541' '542' '543' '544' '545' '546' '547' '548' '549' '550' '551' '552' '553' '554' '555' '556' '557' '558' '559' '560' '561' '562' '563' '564' '565' '566' '567' '568' '569' '570' '571' '572' '573' '574' '575' '576' '577' '578' '579' '580' '581' '582' '583' '584' '585' '586' '587' '588' '589' '590' '591' '592' '593' '594' '595' '596' '597' '598' '599' '600' '601' '602' '603' '604' '605' '606' '607' '608' '609' '610' '611' '612' '613' '614' '615' '616' '617' '618' '619' '620' '621' '622' '623' '624' '625' '626' '627' '628' '629' '630' '631' '632' '633' '634' '635' '636' '637' '638' '639' '640' '641' '642' '643' '644' '645' '646' '647' '648' '649' '650' '651' '652' '653' '654' '655' '656' '657' '658' '659' '660' '661' '662' '663' '664' '665' '666' '667' '668' '669' '670' '671' '672' '673' '674' '675' '676' '677' '678' '679' '680' '681' '682' '683' '684' '685' '686' '687' '688' '689' '690' '691' '692' '693' '694' '695' '696' '697' '698' '699' '700')
# '701' '702' '703' '704' '705' '706' '707' '708' '709' '710' '711' '712' '713' '714' '715' '716' '717' '718' '719' '720' '721' '722' '723' '724' '725' '726' '727' '728' '729' '730' '731' '732' '733' '734' '735' '736' '737' '738' '739' '740' '741' '742' '743' '744' '745' '746' '747' '748' '749' '750' '751' '752' '753' '754' '755' '756' '757' '758' '759' '760' '761' '762' '763' '764' '765' '766' '767' '768' '769' '770' '771' '772' '773' '774' '775' '776' '777' '778' '779' '780' '781' '782' '783' '784' '785' '786' '787' '788' '789' '790' '791' '792' '793' '794' '795' '796' '797' '798' '799')
# '800' '801' '802' '803' '804' '805' '806' '807' '808' '809' '810' '811' '812' '813' '814' '815' '816' '817' '818' '819' '820' '821' '822' '823' '824' '825' '826' '827' '828' '829' '830' '831' '832' '833' '834' '835' '836' '837' '838' '839' '840' '841' '842' '843' '844' '845' '846' '847' '848' '849' '850' '851' '852' '853' '854' '855' '856' '857' '858' '859' '860' '861' '862' '863' '864' '865' '866' '867' '868' '869' '870' '871' '872' '873' '874' '875' '876' '877' '878' '879' '880' '881' '882' '883' '884' '885' '886' '887' '888' '889' '890' '891' '892' '893' '894' '895' '896' '897' '898' '899' '900' '901' '902' '903' '904' '905' '906' '907' '908' '909' '910' '911' '912' '913' '914' '915' '916' '917' '918' '919' '920' '921' '922' '923' '924' '925' '926' '927' '928' '929' '930' '931' '932' '933' '934' '935' '936' '937' '938' '939' '940' '941' '942' '943' '944' '945' '946' '947' '948' '949' '950' '951' '952' '953' '954' '955' '956' '957' '958' '959' '960' '961' '962' '963' '964' '965' '966' '967' '968' '969' '970' '971' '972' '973' '974' '975' '976' '977' '978' '979' '980' '981' '982' '983' '984' '985' '986' '987' '988' '989' '990' '991' '992' '993' '994' '995' '996' '997' '998' '999' '1000' '1001' '1002' '1003' '1004' '1005' '1006' '1007' '1008' '1009' '1010' '1011' '1012' '1013' '1014' '1015' '1016' '1017' '1018' '1019' '1020' '1021' '1022' '1023' '1024' '1025' '1026' '1027' '1028' '1029' '1030' '1031' '1032' '1033' '1034' '1035' '1036' '1037' '1038' '1039' '1040' '1041' '1042' '1043' '1044' '1045' '1046' '1047' '1048' '1049' '1050' '1051' '1052' '1053' '1054' '1055' '1056' '1057' '1058' '1059' '1060' '1061' '1062' '1063' '1064' '1065' '1066' '1067' '1068' '1069' '1070' '1071' '1072' '1073' '1074' '1075' '1076' '1077' '1078' '1079' '1080' '1081' '1082' '1083' '1084' '1085' '1086' '1087' '1088' '1089' '1090' '1091' '1092' '1093' '1094' '1095' '1096' '1097' '1098' '1099' '1100' '1101' '1102' '1103' '1104' '1105' '1106' '1107' '1108' '1109' '1110' '1111' '1112' '1113' '1114' '1115' '1116' '1117' '1118' '1119' '1120' '1121' '1122' '1123' '1124' '1125' '1126' '1127' '1128' '1129' '1130' '1131' '1132' '1133' '1134' '1135' '1136' '1137' '1138' '1139' '1140' '1141' '1142' '1143' '1144' '1145' '1146' '1147' '1148' '1149' '1150' '1151' '1152' '1153' '1154' '1155' '1156' '1157' '1158' '1159' '1160' '1161' '1162' '1163' '1164' '1165' '1166' '1167' '1168' '1169' '1170' '1171' '1172' '1173' '1174' '1175' '1176' '1177' '1178' '1179' '1180' '1181' '1182' '1183' '1184' '1185' '1186' '1187' '1188' '1189' '1190' '1191' '1192' '1193' '1194' '1195' '1196' '1197' '1198' '1199' '1200' '1201' '1202' '1203' '1204' '1205' '1206' '1207' '1208' '1209' '1210' '1211' '1212' '1213' '1214' '1215' '1216' '1217' '1218' '1219' '1220' '1221' '1222' '1223' '1224' '1225' '1226' '1227' '1228' '1229' '1230' '1231' '1232' '1233' '1234' '1235' '1236' '1237' '1238' '1239' '1240' '1241' '1242' '1243' '1244' '1245' '1246' '1247' '1248' '1249' '1250' '1251' '1252' '1253' '1254' '1255' '1256' '1257' '1258' '1259' '1260' '1261' '1262' '1263' '1264' '1265' '1266' '1267' '1268' '1269' '1270' '1271' '1272' '1273' '1274' '1275' '1276' '1277' '1278' '1279' '1280' '1281' '1282' '1283' '1284' '1285' '1286' '1287' '1288' '1289' '1290' '1291' '1292' '1293' '1294' '1295' '1296' '1297' '1298' '1299' '1300' '1301' '1302' '1303' '1304' '1305' '1306' '1307' '1308' '1309' '1310' '1311' '1312' '1313' '1314' '1315' '1316' '1317' '1318' '1319' '1320' '1321' '1322' '1323' '1324' '1325' '1326' '1327' '1328' '1329' '1330' '1331' '1332' '1333' '1334' '1335' '1336' '1337' '1338' '1339' '1340' '1341' '1342' '1343' '1344' '1345' '1346' '1347' '1348' '1349' '1350' '1351' '1352' '1353' '1354' '1355' '1356' '1357' '1358' '1359' '1360' '1361' '1362' '1363' '1364' '1365' '1366' '1367' '1368' '1369' '1370' '1371' '1372' '1373' '1374' '1375' '1376' '1377' '1378' '1379' '1380' '1381' '1382' '1383' '1384' '1385' '1386' '1387' '1388' '1389' '1390' '1391' '1392' '1393' '1394' '1395' '1396' '1397' '1398' '1399' '1400')

set species = ('R0001' 'R0002' 'R0003' 'R0004' 'R0005' 'R0006' 'R0007' 'R0008' 'R0009' 'R0010' 'R0011' 'R0012' 'R0013' 'R0014' 'R0015' 'R0016' 'R0017' 'R0018' 'R0019' 'R0020' 'R0021' 'R0022' 'R0023' 'R0024' 'R0025' 'R0026' 'R0027' 'R0028' 'R0029' 'R0030' 'R0031' 'R0032' 'R0033' 'R0034' 'R0035' 'R0036' 'R0037' 'R0038' 'R0039' 'R0040' 'R0041' 'R0042' 'R0043' 'R0044' 'R0045' 'R0046' 'R0047' 'R0048' 'R0049' 'R0050' 'R0051' 'R0052' 'R0053' 'R0054' 'R0055' 'R0056' 'R0057' 'R0058' 'R0059' 'R0060' 'R0061' 'R0062' 'R0063' 'R0064' 'R0065' 'R0066' 'R0067' 'R0068' 'R0069' 'R0070' 'R0071' 'R0072' 'R0073' 'R0074' 'R0075' 'R0076' 'R0077' 'R0078' 'R0079' 'R0080' 'R0081' 'R0082' 'R0083' 'R0084' 'R0085' 'R0086' 'R0087' 'R0088' 'R0089' 'R0090' 'R0091' 'R0092' 'R0093' 'R0094' 'R0095' 'R0096' 'R0097' 'R0098' 'R0099' 'R0100' 'R0101' 'R0102' 'R0103' 'R0104' 'R0105' 'R0106' 'R0107' 'R0108' 'R0109' 'R0110' 'R0111' 'R0112' 'R0113' 'R0114' 'R0115' 'R0116' 'R0117' 'R0118' 'R0119' 'R0120' 'R0121' 'R0122' 'R0123' 'R0124' 'R0125' 'R0126' 'R0127' 'R0128' 'R0129' 'R0130' 'R0131' 'R0132' 'R0133' 'R0134' 'R0135' 'R0136' 'R0137' 'R0138' 'R0139' 'R0140' 'R0141' 'R0142' 'R0143' 'R0144' 'R0145' 'R0146' 'R0147' 'R0148' 'R0149' 'R0150' 'R0151' 'R0152' 'R0153' 'R0154' 'R0155' 'R0156' 'R0157' 'R0158' 'R0159' 'R0160' 'R0161' 'R0162' 'R0163' 'R0164' 'R0165' 'R0166' 'R0167' 'R0168' 'R0169' 'R0170' 'R0171' 'R0172' 'R0173' 'R0174' 'R0175' 'R0176' 'R0177' 'R0178' 'R0179' 'R0180' 'R0181' 'R0182' 'R0183' 'R0184' 'R0185' 'R0186' 'R0187' 'R0188' 'R0189' 'R0190' 'R0191' 'R0192' 'R0193' 'R0194' 'R0195' 'R0196' 'R0197' 'R0198' 'R0199' 'R0200' 'R0201' 'R0202' 'R0203' 'R0204' 'R0205' 'R0206' 'R0207' 'R0208' 'R0209' 'R0210' 'R0211' 'R0212' 'R0213' 'R0214' 'R0215' 'R0216' 'R0217' 'R0218' 'R0219' 'R0220' 'R0221' 'R0222' 'R0223' 'R0224' 'R0225' 'R0226' 'R0227' 'R0228' 'R0229' 'R0230' 'R0231' 'R0232' 'R0233' 'R0234' 'R0235' 'R0236' 'R0237' 'R0238' 'R0239' 'R0240' 'R0241' 'R0242' 'R0243' 'R0244' 'R0245' 'R0246' 'R0247' 'R0248' 'R0249' 'R0250' 'R0251' 'R0252' 'R0253' 'R0254' 'R0255' 'R0256' 'R0257' 'R0258' 'R0259' 'R0260' 'R0261' 'R0262' 'R0263' 'R0264' 'R0265' 'R0266' 'R0267' 'R0268' 'R0269' 'R0270' 'R0271' 'R0272' 'R0273' 'R0274' 'R0275' 'R0276' 'R0277' 'R0278' 'R0279' 'R0280' 'R0281' 'R0282' 'R0283' 'R0284' 'R0285' 'R0286' 'R0287' 'R0288' 'R0289' 'R0290' 'R0291' 'R0292' 'R0293' 'R0294' 'R0295' 'R0296' 'R0297' 'R0298' 'R0299' 'R0300' 'R0301' 'R0302' 'R0303' 'R0304' 'R0305' 'R0306' 'R0307' 'R0308' 'R0309' 'R0310' 'R0311' 'R0312' 'R0313' 'R0314' 'R0315' 'R0316' 'R0317' 'R0318' 'R0319' 'R0320' 'R0321' 'R0322' 'R0323' 'R0324' 'R0325' 'R0326' 'R0327' 'R0328' 'R0329' 'R0330' 'R0331' 'R0332' 'R0333' 'R0334' 'R0335' 'R0336' 'R0337' 'R0338' 'R0339' 'R0340' 'R0341' 'R0342' 'R0343' 'R0344' 'R0345' 'R0346' 'R0347' 'R0348' 'R0349' 'R0350' 'R0351' 'R0352' 'R0353' 'R0354' 'R0355' 'R0356' 'R0357' 'R0358' 'R0359' 'R0360' 'R0361' 'R0362' 'R0363' 'R0364' 'R0365' 'R0366' 'R0367' 'R0368' 'R0369' 'R0370' 'R0371' 'R0372' 'R0373' 'R0374' 'R0375' 'R0376' 'R0377' 'R0378' 'R0379' 'R0380' 'R0381' 'R0382' 'R0383' 'R0384' 'R0385' 'R0386' 'R0387' 'R0388' 'R0389' 'R0390' 'R0391' 'R0392' 'R0393' 'R0394' 'R0395' 'R0396' 'R0397' 'R0398' 'R0399' 'R0400' 'R0401' 'R0402' 'R0403' 'R0404' 'R0405' 'R0406' 'R0407' 'R0408' 'R0409' 'R0410' 'R0411' 'R0412' 'R0413' 'R0414' 'R0415' 'R0416' 'R0417' 'R0418' 'R0419' 'R0420' 'R0421' 'R0422' 'R0423' 'R0424' 'R0425' 'R0426' 'R0427' 'R0428' 'R0429' 'R0430' 'R0431' 'R0432' 'R0433' 'R0434' 'R0435' 'R0436' 'R0437' 'R0438' 'R0439' 'R0440')
#'R0441' 'R0442' 'R0443' 'R0444' 'R0445' 'R0446' 'R0447' 'R0448' 'R0449' 'R0450' 'R0451' 'R0452' 'R0453' 'R0454' 'R0455' 'R0456' 'R0457' 'R0458' 'R0459' 'R0460' 'R0461' 'R0462' 'R0463' 'R0464' 'R0465' 'R0466' 'R0467' 'R0468' 'R0469' 'R0470' 'R0471' 'R0472' 'R0473' 'R0474' 'R0475' 'R0476' 'R0477' 'R0478' 'R0479' 'R0480' 'R0481' 'R0482' 'R0483' 'R0484' 'R0485' 'R0486' 'R0487' 'R0488' 'R0489' 'R0490' 'R0491' 'R0492' 'R0493' 'R0494' 'R0495' 'R0496' 'R0497' 'R0498' 'R0499' 'R0500' 'R0501' 'R0502' 'R0503' 'R0504' 'R0505' 'R0506' 'R0507' 'R0508' 'R0509' 'R0510' 'R0511' 'R0512' 'R0513' 'R0514' 'R0515' 'R0516' 'R0517' 'R0518' 'R0519' 'R0520' 'R0521' 'R0522' 'R0523' 'R0524' 'R0525' 'R0526' 'R0527' 'R0528' 'R0529' 'R0530' 'R0531' 'R0532' 'R0533' 'R0534' 'R0535' 'R0536' 'R0537' 'R0538' 'R0539' 'R0540' 'R0541' 'R0542' 'R0543' 'R0544' 'R0545' 'R0546' 'R0547' 'R0548' 'R0549' 'R0550' 'R0551' 'R0552' 'R0553' 'R0554' 'R0555' 'R0556' 'R0557' 'R0558' 'R0559' 'R0560' 'R0561' 'R0562' 'R0563' 'R0564' 'R0565' 'R0566' 'R0567' 'R0568' 'R0569' 'R0570' 'R0571' 'R0572' 'R0573' 'R0574' 'R0575' 'R0576' 'R0577' 'R0578' 'R0579' 'R0580' 'R0581' 'R0582' 'R0583' 'R0584' 'R0585' 'R0586' 'R0587' 'R0588' 'R0589' 'R0590' 'R0591' 'R0592' 'R0593' 'R0594' 'R0595' 'R0596' 'R0597' 'R0598' 'R0599' 'R0600' 'R0601' 'R0602' 'R0603' 'R0604' 'R0605' 'R0606' 'R0607' 'R0608' 'R0609' 'R0610' 'R0611' 'R0612' 'R0613' 'R0614' 'R0615' 'R0616' 'R0617' 'R0618' 'R0619' 'R0620' 'R0621' 'R0622' 'R0623' 'R0624' 'R0625' 'R0626' 'R0627' 'R0628' 'R0629' 'R0630' 'R0631' 'R0632' 'R0633' 'R0634' 'R0635' 'R0636' 'R0637' 'R0638' 'R0639' 'R0640' 'R0641' 'R0642' 'R0643' 'R0644' 'R0645' 'R0646' 'R0647' 'R0648' 'R0649' 'R0650' 'R0651' 'R0652' 'R0653' 'R0654' 'R0655' 'R0656' 'R0657' 'R0658' 'R0659' 'R0660' 'R0661' 'R0662' 'R0663' 'R0664' 'R0665' 'R0666' 'R0667' 'R0668' 'R0669' 'R0670' 'R0671' 'R0672' 'R0673' 'R0674' 'R0675' 'R0676' 'R0677' 'R0678' 'R0679' 'R0680' 'R0681' 'R0682' 'R0683' 'R0684' 'R0685' 'R0686' 'R0687' 'R0688' 'R0689' 'R0690' 'R0691' 'R0692' 'R0693' 'R0694' 'R0695' 'R0696' 'R0697' 'R0698' 'R0699' 'R0700')
# 'R0701' 'R0702' 'R0703' 'R0704' 'R0705' 'R0706' 'R0707' 'R0708' 'R0709' 'R0710' 'R0711' 'R0712' 'R0713' 'R0714' 'R0715' 'R0716' 'R0717' 'R0718' 'R0719' 'R0720' 'R0721' 'R0722' 'R0723' 'R0724' 'R0725' 'R0726' 'R0727' 'R0728' 'R0729' 'R0730' 'R0731' 'R0732' 'R0733' 'R0734' 'R0735' 'R0736' 'R0737' 'R0738' 'R0739' 'R0740' 'R0741' 'R0742' 'R0743' 'R0744' 'R0745' 'R0746' 'R0747' 'R0748' 'R0749' 'R0750' 'R0751' 'R0752' 'R0753' 'R0754' 'R0755' 'R0756' 'R0757' 'R0758' 'R0759' 'R0760' 'R0761' 'R0762' 'R0763' 'R0764' 'R0765' 'R0766' 'R0767' 'R0768' 'R0769' 'R0770' 'R0771' 'R0772' 'R0773' 'R0774' 'R0775' 'R0776' 'R0777' 'R0778' 'R0779' 'R0780' 'R0781' 'R0782' 'R0783' 'R0784' 'R0785' 'R0786' 'R0787' 'R0788' 'R0789' 'R0790' 'R0791' 'R0792' 'R0793' 'R0794' 'R0795' 'R0796' 'R0797' 'R0798' 'R0799')
# 'R0800' 'R0801' 'R0802' 'R0803' 'R0804' 'R0805' 'R0806' 'R0807' 'R0808' 'R0809' 'R0810' 'R0811' 'R0812' 'R0813' 'R0814' 'R0815' 'R0816' 'R0817' 'R0818' 'R0819' 'R0820' 'R0821' 'R0822' 'R0823' 'R0824' 'R0825' 'R0826' 'R0827' 'R0828' 'R0829' 'R0830' 'R0831' 'R0832' 'R0833' 'R0834' 'R0835' 'R0836' 'R0837' 'R0838' 'R0839' 'R0840' 'R0841' 'R0842' 'R0843' 'R0844' 'R0845' 'R0846' 'R0847' 'R0848' 'R0849' 'R0850' 'R0851' 'R0852' 'R0853' 'R0854' 'R0855' 'R0856' 'R0857' 'R0858' 'R0859' 'R0860' 'R0861' 'R0862' 'R0863' 'R0864' 'R0865' 'R0866' 'R0867' 'R0868' 'R0869' 'R0870' 'R0871' 'R0872' 'R0873' 'R0874' 'R0875' 'R0876' 'R0877' 'R0878' 'R0879' 'R0880' 'R0881' 'R0882' 'R0883' 'R0884' 'R0885' 'R0886' 'R0887' 'R0888' 'R0889' 'R0890' 'R0891' 'R0892' 'R0893' 'R0894' 'R0895' 'R0896' 'R0897' 'R0898' 'R0899' 'R0900' 'R0901' 'R0902' 'R0903' 'R0904' 'R0905' 'R0906' 'R0907' 'R0908' 'R0909' 'R0910' 'R0911' 'R0912' 'R0913' 'R0914' 'R0915' 'R0916' 'R0917' 'R0918' 'R0919' 'R0920' 'R0921' 'R0922' 'R0923' 'R0924' 'R0925' 'R0926' 'R0927' 'R0928' 'R0929' 'R0930' 'R0931' 'R0932' 'R0933' 'R0934' 'R0935' 'R0936' 'R0937' 'R0938' 'R0939' 'R0940' 'R0941' 'R0942' 'R0943' 'R0944' 'R0945' 'R0946' 'R0947' 'R0948' 'R0949' 'R0950' 'R0951' 'R0952' 'R0953' 'R0954' 'R0955' 'R0956' 'R0957' 'R0958' 'R0959' 'R0960' 'R0961' 'R0962' 'R0963' 'R0964' 'R0965' 'R0966' 'R0967' 'R0968' 'R0969' 'R0970' 'R0971' 'R0972' 'R0973' 'R0974' 'R0975' 'R0976' 'R0977' 'R0978' 'R0979' 'R0980' 'R0981' 'R0982' 'R0983' 'R0984' 'R0985' 'R0986' 'R0987' 'R0988' 'R0989' 'R0990' 'R0991' 'R0992' 'R0993' 'R0994' 'R0995' 'R0996' 'R0997' 'R0998' 'R0999' 'R1000' 'R1001' 'R1002' 'R1003' 'R1004' 'R1005' 'R1006' 'R1007' 'R1008' 'R1009' 'R1010' 'R1011' 'R1012' 'R1013' 'R1014' 'R1015' 'R1016' 'R1017' 'R1018' 'R1019' 'R1020' 'R1021' 'R1022' 'R1023' 'R1024' 'R1025' 'R1026' 'R1027' 'R1028' 'R1029' 'R1030' 'R1031' 'R1032' 'R1033' 'R1034' 'R1035' 'R1036' 'R1037' 'R1038' 'R1039' 'R1040' 'R1041' 'R1042' 'R1043' 'R1044' 'R1045' 'R1046' 'R1047' 'R1048' 'R1049' 'R1050' 'R1051' 'R1052' 'R1053' 'R1054' 'R1055' 'R1056' 'R1057' 'R1058' 'R1059' 'R1060' 'R1061' 'R1062' 'R1063' 'R1064' 'R1065' 'R1066' 'R1067' 'R1068' 'R1069' 'R1070' 'R1071' 'R1072' 'R1073' 'R1074' 'R1075' 'R1076' 'R1077' 'R1078' 'R1079' 'R1080' 'R1081' 'R1082' 'R1083' 'R1084' 'R1085' 'R1086' 'R1087' 'R1088' 'R1089' 'R1090' 'R1091' 'R1092' 'R1093' 'R1094' 'R1095' 'R1096' 'R1097' 'R1098' 'R1099' 'R1100' 'R1101' 'R1102' 'R1103' 'R1104' 'R1105' 'R1106' 'R1107' 'R1108' 'R1109' 'R1110' 'R1111' 'R1112' 'R1113' 'R1114' 'R1115' 'R1116' 'R1117' 'R1118' 'R1119' 'R1120' 'R1121' 'R1122' 'R1123' 'R1124' 'R1125' 'R1126' 'R1127' 'R1128' 'R1129' 'R1130' 'R1131' 'R1132' 'R1133' 'R1134' 'R1135' 'R1136' 'R1137' 'R1138' 'R1139' 'R1140' 'R1141' 'R1142' 'R1143' 'R1144' 'R1145' 'R1146' 'R1147' 'R1148' 'R1149' 'R1150' 'R1151' 'R1152' 'R1153' 'R1154' 'R1155' 'R1156' 'R1157' 'R1158' 'R1159' 'R1160' 'R1161' 'R1162' 'R1163' 'R1164' 'R1165' 'R1166' 'R1167' 'R1168' 'R1169' 'R1170' 'R1171' 'R1172' 'R1173' 'R1174' 'R1175' 'R1176' 'R1177' 'R1178' 'R1179' 'R1180' 'R1181' 'R1182' 'R1183' 'R1184' 'R1185' 'R1186' 'R1187' 'R1188' 'R1189' 'R1190' 'R1191' 'R1192' 'R1193' 'R1194' 'R1195' 'R1196' 'R1197' 'R1198' 'R1199' 'R1200' 'R1201' 'R1202' 'R1203' 'R1204' 'R1205' 'R1206' 'R1207' 'R1208' 'R1209' 'R1210' 'R1211' 'R1212' 'R1213' 'R1214' 'R1215' 'R1216' 'R1217' 'R1218' 'R1219' 'R1220' 'R1221' 'R1222' 'R1223' 'R1224' 'R1225' 'R1226' 'R1227' 'R1228' 'R1229' 'R1230' 'R1231' 'R1232' 'R1233' 'R1234' 'R1235' 'R1236' 'R1237' 'R1238' 'R1239' 'R1240' 'R1241' 'R1242' 'R1243' 'R1244' 'R1245' 'R1246' 'R1247' 'R1248' 'R1249' 'R1250' 'R1251' 'R1252' 'R1253' 'R1254' 'R1255' 'R1256' 'R1257' 'R1258' 'R1259' 'R1260' 'R1261' 'R1262' 'R1263' 'R1264' 'R1265' 'R1266' 'R1267' 'R1268' 'R1269' 'R1270' 'R1271' 'R1272' 'R1273' 'R1274' 'R1275' 'R1276' 'R1277' 'R1278' 'R1279' 'R1280' 'R1281' 'R1282' 'R1283' 'R1284' 'R1285' 'R1286' 'R1287' 'R1288' 'R1289' 'R1290' 'R1291' 'R1292' 'R1293' 'R1294' 'R1295' 'R1296' 'R1297' 'R1298' 'R1299' 'R1300' 'R1301' 'R1302' 'R1303' 'R1304' 'R1305' 'R1306' 'R1307' 'R1308' 'R1309' 'R1310' 'R1311' 'R1312' 'R1313' 'R1314' 'R1315' 'R1316' 'R1317' 'R1318' 'R1319' 'R1320' 'R1321' 'R1322' 'R1323' 'R1324' 'R1325' 'R1326' 'R1327' 'R1328' 'R1329' 'R1330' 'R1331' 'R1332' 'R1333' 'R1334' 'R1335' 'R1336' 'R1337' 'R1338' 'R1339' 'R1340' 'R1341' 'R1342' 'R1343' 'R1344' 'R1345' 'R1346' 'R1347' 'R1348' 'R1349' 'R1350' 'R1351' 'R1352' 'R1353' 'R1354' 'R1355' 'R1356' 'R1357' 'R1358' 'R1359' 'R1360' 'R1361' 'R1362' 'R1363' 'R1364' 'R1365' 'R1366' 'R1367' 'R1368' 'R1369' 'R1370' 'R1371' 'R1372' 'R1373' 'R1374' 'R1375' 'R1376' 'R1377' 'R1378' 'R1379' 'R1380' 'R1381' 'R1382' 'R1383' 'R1384' 'R1385' 'R1386' 'R1387' 'R1388' 'R1389' 'R1390' 'R1391' 'R1392' 'R1393' 'R1394' 'R1395' 'R1396' 'R1397' 'R1398' 'R1399' 'R1400')

# just exchange
# set n_spec = ('1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' '23' '24' '25' '26' '27')
# set species = ('Cl2ex1' 'HOClex1' 'HClex1' 'BrClex1' 'HOBrex1' 'HBrex1' 'Br2ex1' 'Cl2ex2' 'HOClex2' 'HClex2' 'BrClex2' 'HOBrex2' 'HBrex2' 'Br2ex2' 'Cl2ex3' 'HOClex3' 'HClex3' 'BrClex3' 'HOBrex3' 'HBrex3' 'Br2ex3' 'Clnet1' 'Brnet1' 'Clnet2' 'Brnet2' 'Clnet3' 'Brnet3')

# ===========================================================================
# HANDS OFF FROM BELOW!! NO USER SERVICABLE PARTS INSIDE!! DEVELOPPERS ONLY!!
# ===========================================================================

# make sure that each plotprog file is saved and none overwritten to be able to
# reproduce runs; this plotprog file name is put into plot

# numbers of plots per page
set pl_n = ("7")

# define "overflow" value for viewports
set pl_np = $pl_n
@ pl_np+=1

# define header
set n_hl=("1" "2" "3" "4" "5" "6" "7" "8")
set hl1=("! Description: Ferret plotprogram for Mistra netCDf output: reaction rates\n")
set hl2=("! input unit is mol m-3 s-1\n")
set hl3=(" \n")
set hl4=("! reset everything\n")
set hl5=("\cancel mode verify ! like "unset echo" under unix\n")
set hl6=("cancel data/all\n")
set hl7=("cancel variable/all\n")
set hl8=("cancel symbol/all\n")

# write header
echo $hl1 $hl2 $hl3 $hl4 $hl6 $hl7 $hl8  >> pgtmp0

# write special definitions
echo "let Cl2ex1=r0458-r0459"  >> pgtmp0
echo "let Cl2ex2=r0739-r0740"  >> pgtmp0
echo "let Cl2ex3=r1056-r1057"  >> pgtmp0
echo "let HOClex1=r0456-r0457"  >> pgtmp0
echo "let HOClex2=r0737-r0738"  >> pgtmp0
echo "let HOClex3=r1054-r1055"  >> pgtmp0
echo "let HClex1=r0454-r0455"  >> pgtmp0
echo "let HClex2=r0735-r0736"  >> pgtmp0
echo "let HClex3=r1052-r1053"  >> pgtmp0
echo "let BrClex1=r0466-r0467"  >> pgtmp0
echo "let BrClex2=r0747-r0748"  >> pgtmp0
echo "let BrClex3=r1064-r1065"  >> pgtmp0
echo "let HOBrex1=r0462-r0463"  >> pgtmp0
echo "let HOBrex2=r0743-r0744"  >> pgtmp0
echo "let HOBrex3=r1060-r1061"  >> pgtmp0
echo "let HBrex1=r0460-r0461"  >> pgtmp0
echo "let HBrex2=r0741-r0742"  >> pgtmp0
echo "let HBrex3=r1058-r1059"  >> pgtmp0
echo "let Br2ex1=r0464-r0465"  >> pgtmp0
echo "let Br2ex2=r0745-r0746"  >> pgtmp0
echo "let Br2ex3=r1062-r1063"  >> pgtmp0

echo "let Clnet1=2*Cl2ex1+HOClex1+HClex1+BrClex1-r0208+r0210+r0211+r0212-r0216"   >> pgtmp0
echo "let Clnet2=2*Cl2ex2+HOClex2+HClex2+BrClex2-r0226+r0228+r0229+r0229-r0234"   >> pgtmp0
echo "let Clnet3=2*Cl2ex3+HOClex3+HClex3+BrClex3-r0806+r0808+r0809+r0810-r0814"   >> pgtmp0
echo "let Brnet1=2*Br2ex1+HOBrex1+HBrex1+BrClex1-r0209+r0213+r0214+r0215-r0217-r0218"   >> pgtmp0
echo "let Brnet2=2*Br2ex2+HOBrex2+HBrex2+BrClex2-r0227+r0231+r0232+r0233-r0235-r0236"   >> pgtmp0
echo "let Brnet3=2*Br2ex3+HOBrex3+HBrex3+BrClex3-r0807+r0811+r0812+r0813-r0815-r0816"   >> pgtmp0


# open files
foreach i ($m_runs)
    echo "use" ' "'  "$runs[$i]"'rxnrate.nc"' '; use "'  "$runs[$i]"'meteo.nc"' >> pgtmp0
end

# define timeaxis
# mtime as in netCDF is "local time", get "time since model start"
echo "let mmtime=mtime-mtime[l=1]" >> pgtmp0

# define viewports and page style
# 3x4
echo "define viewport /xlim=0.00,0.49 /ylim=0.65,1.00 VP1"  >> pgtmp0
echo "define viewport /xlim=0.00,0.49 /ylim=0.45,0.80 VP2"  >> pgtmp0
echo "define viewport /xlim=0.00,0.49 /ylim=0.25,0.60 VP3"  >> pgtmp0
echo "define viewport /xlim=0.00,0.49 /ylim=0.05,0.40 VP4"  >> pgtmp0
echo "define viewport /xlim=0.51,1.00 /ylim=0.65,1.00 VP5"  >> pgtmp0
echo "define viewport /xlim=0.51,1.00 /ylim=0.45,0.80 VP6"  >> pgtmp0
echo "define viewport /xlim=0.51,1.00 /ylim=0.25,0.60 VP7"  >> pgtmp0
echo "define viewport /xlim=0.51,1.00 /ylim=0.05,0.40 VP8"  >> pgtmp0

# output in metafile format
echo "SET MODE METAFILE:"$pl_name.1".plt" >> pgtmp0

# plot command is being composed of the different single steps:

# base plot command
if ($pl_ty == 1) then
# 1 - evolution w/ time;
 set plcmd = ("plot /vs /line /nolabel /k=")
 set plcmdo = ("plot /vs /line /nolabel /overlay /k=")
 set pll1 = ("k=")
 set lim1 = ("/vlimits=0:")
 set lim2 = (",l=@max")
 set xax = ("mmtimedd2 ,")
 set yax = (" ")
endif
if ($pl_ty == 2) then
# 2 - vertical profile
# set plcmd = ("plot /nolabel /k=$minheight:$maxheight /l=")
 set plcmd = ("plot /vs /line /nolabel /k=2:$maxheight /l=")
# set plcmdo = ("plot /nolabel /overlay /k=$minheight:$maxheight /l=")
 set plcmdo = ("plot /vs /line /nolabel /overlay /k=2:$maxheight /l=")
 set pll1 = ("l=")
 set lim1 = ("/hlimits=0:")
# set lim2 = (",k=$minheight:$maxheight@max")
 set lim2 = (",k=2:$maxheight@max")
 set xax = (" ")
 set yax = (",etadd2")
endif
if ($pl_ty == 3) then
# 3 - contour plot
# set plcmd = ("shade /nolabel /k=$minheight:")
 set plcmd = ("shade /nolabel /k=2:")
# set plcmdo = ("shade /nolabel /overlay /k=$minheight:")
 set plcmdo = ("shade /nolabel /overlay /k=2:")
 set pll1 = ("k=")
 set lim1 = (" ")
 set lim2 = (",l=@max")
 set pld=$maxheight
 set n_pl = ("1")
 set m_pl = ("1")
 set n_runs = ("1")
 set m_runs = ("1")
 set xax = (" ")
 set yax = (" ")
endif
# 4 - box model
if ($pl_ty == 4) then
# like evolution w/ time, but /k=2
 set plcmd = ("plot /nolabel /k=")
 set plcmdo = ("plot /nolabel /overlay /k=")
 set pll1 = ("k=")
 set lim1 = ("/vlimits=0:")
 set lim2 = (",l=@max")
 set pld=("1")
 set n_pl = ("1")
 set m_pl = ("1")
 set xax = (" ")
 set yax = (" ")
endif
# /title='"'$species[$k]'"'

# unit conversion: 1 - mol m-3 s-1; 2 - mol mol-1 s-1; 3 - molec cm-3 s-1
# unchanged
if ($pl_un == 1) then
    set unit1a=1.
    set unit1b=1.
    set unit1c=1.
    set unit3a=1.
    set unit3b=1.
    set unit3c=1.
    set unit5a=1.
    set unit5b=1.
    set unit5c=1.
    set unit7a=1.
    set unit7b=1.
    set unit7c=1.
endif
# mol m-3 --> mol mol-1: 1./(rho / 29.e-3)
if ($pl_un == 2) then
    set unit1a=34.4827
    set unit1b="rho[de2,$pll1$pld[1]$lim2]"
    set unit1c="rhodd2"
    set unit3a=34.4827
    set unit3b="rho[de4,$pll1$pld[1]$lim2]"
    set unit3c="rhodd4"
    set unit5a=34.4827
    set unit5b="rho[de6,$pll1$pld[1]$lim2]"
    set unit5c="rhodd6"
    set unit7a=34.4827
    set unit7b="rho[de6,$pll1$pld[1]$lim2]"
    set unit7c="rhodd8"
endif
# mol m-3 --> molec cm-3: 1./(1.e6 / 6.022e23)
if ($pl_un == 3) then
    set unit1a=1.6606e-18
    set unit1b=1.
    set unit1c=1.
    set unit3a=1.6606e-18
    set unit3b=1.
    set unit3c=1.
    set unit5a=1.6606e-18
    set unit5b=1.
    set unit5c=1.
    set unit7a=1.6606e-18
    set unit7b=1.
    set unit7c=1.
endif

# run info always on view port1 (VP1)
set countVP = 1
echo "set viewport VP$countVP" >> pgtmp0
# define linestyle master
set lstym=("" "/dash=(.1,.1,.1,.1)" "/dash=(.01,.1,.01,.1)" "/dash=(.3,.1,.3,.1)")
set lsty=("" "" "" "" "" "" "" "" "" "")
# define labels
set line1 =(" ")
set line2 =(" ")
set line3 =(" ")
set line4 =(" ")
set line5 =(" ")
set line6 =(" ")
set line7 =(" ")
set line8 =(" ")
set line9 =(" ")

# shorten directory names, max 32. chars
# n_runs - number of runs
# runs - complete directory names of runs
# TODO! maybe via awk and substr command - output inextra file, but how to read in from that file into var?

# one line explanation of run
if ($n_runs == 1) then
    set countl = 0
    foreach kk ($m_pl)
       @ countl+=1
       set line$countl =("$runs[1] $pld[$kk]")
       set lsty[$countl] = $lstym[$kk]
    end
endif
if ($n_runs == 2) then
    set countl = 0
    foreach kk ($m_pl)
       @ countl+=1
       set line$countl =("$runs[1] $pld[$kk]")
       set lsty[$countl] = $lstym[$kk]
       @ countl+=1
       set line$countl =("$runs[2] $pld[$kk]")
       set lsty[$countl] = $lstym[$kk]
    end
endif
if ($n_runs == 3) then
    set countl = 0
    foreach kk ($m_pl)
       @ countl+=1
       set line$countl =("$runs[1] $pld[$kk]")
       set lsty[$countl] = $lstym[$kk]
       @ countl+=1
       set line$countl =("$runs[2] $pld[$kk]")
       set lsty[$countl] = $lstym[$kk]
       @ countl+=1
       set line$countl =("$runs[3] $pld[$kk]")
       set lsty[$countl] = $lstym[$kk]
    end
endif
if ($n_runs == 4) then
    set countl = 0
    foreach kk ($m_pl)
       @ countl+=1
       set line$countl =("$runs[1] $pld[$kk]")
       set lsty[$countl] = $lstym[$kk]
       @ countl+=1
       set line$countl =("$runs[2] $pld[$kk]")
       set lsty[$countl] = $lstym[$kk]
       @ countl+=1
       set line$countl =("$runs[3] $pld[$kk]")
       set lsty[$countl] = $lstym[$kk]
       @ countl+=1
       set line$countl =("$runs[4] $pld[$kk]")
       set lsty[$countl] = $lstym[$kk]
    end
endif

# put unit on plot
# set unit_c=(" ")
# if ($pl_un == 1) set unit_c=("unit: mol m-3")
# if ($pl_un == 2) set unit_c=("mol mol1")
# if ($pl_un == 3) set unit_c=("unit: molec cm-3")
set unit_c=("mol/(m3_s)" "mol/(mol_s)" "molec/(cm3_s)")
set line9 = $unit_c[$pl_un]

# empty plot
echo "plot /i=0:100 /hlimits=0:100 /vlimits=0:100 /noaxis /nolabel 0" >> pgtmp0
# plot explanation
echo "label 10,110,  -1, 0, .25 @P1Gas phase" >> pgtmp0
echo "label 0, 90, -1, 0, .15 @P1$line1" >> pgtmp0
echo "plot /overlay /vs /nolabel /line=1 $lsty[1] {80,99,99}, {95,95,95}" >> pgtmp0
echo "label 0, 80, -1, 0, .15 @P1$line2" >> pgtmp0
if ($countl >= 2) echo "plot /overlay /vs /nolabel /line=2 $lsty[2] {80,99,99}, {85,85,85}" >> pgtmp0
echo "label 0, 70, -1, 0, .15 @P1$line3" >> pgtmp0
if ($countl >= 3) echo "plot /overlay /vs /nolabel /line=3 $lsty[3] {80,99,99}, {75,75,75}" >> pgtmp0
echo "label 0, 60, -1, 0, .15 @P1$line4" >> pgtmp0
if ($countl >= 4) echo "plot /overlay /vs /nolabel /line=4 $lsty[4] {80,99,99}, {65,65,65}" >> pgtmp0
echo "label 0, 50, -1, 0, .15 @P1$line5" >> pgtmp0
if ($countl >= 5) echo "plot /overlay /vs /nolabel /line=5 $lsty[5] {80,99,99}, {55,55,55}" >> pgtmp0
echo "label 0, 40, -1, 0, .15 @P1$line6" >> pgtmp0
if ($countl >= 6) echo "plot /overlay /vs /nolabel /line=6 $lsty[6] {80,99,99}, {45,45,45}" >> pgtmp0
echo "label 0, 30, -1, 0, .15 @P1$line7" >> pgtmp0
if ($countl >= 7) echo "plot /overlay /vs /nolabel /line=7 $lsty[7] {80,99,99}, {35,35,35}" >> pgtmp0
echo "label 0, 20, -1, 0, .15 @P1$line8" >> pgtmp0
if ($countl >= 8) echo "plot /overlay /vs /nolabel /line=8 $lsty[8] {80,99,99}, {25,25,25}" >> pgtmp0
echo "label 0, 10, -1, 0, .15 @P1$line9" >> pgtmp0
if ($countl >= 9) echo "plot /overlay /vs /nolabel /line=9 $lsty[9] {80,99,99}, {15,15,15}" >> pgtmp0

# counter for pages
set countPG = 1

# plot loop
foreach k ($n_spec)
# number of runs to be overplotted (= n_runs)
    @ countVP+=1
    echo "set viewport VP$countVP" >> pgtmp0
# counter, set viewport w/ counter: set viewport VP$count
# d1, d3, d5, .. have to be replaced with sed in final ferret file to give [d=1] etc

# echo "after 8"
# echo $k
# echo $n_spec

# 1 run
    if ($n_runs == 1) then
#       find max to correctly scale axes
#        echo "let vmax=0." >> pgtmp0
#       foreach kk($m_pl)
#           echo "let vmax1=1.05*$species[$k][de1,$pll1$pld[$kk]$lim2]/$unit1a/$unit1b" >> pgtmp0
#           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
#       end
#       echo "if "'`'"vmax eq 0."'`'" then let vmax=1 endif" >> pgtmp0
##       plot
#       if ($pl_ty != 3) then
#           echo $plcmd$pld[1] $lstym[1] $lim1'`'"1.05*vmax"'`'  $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
#       else
#           echo $plcmd$pld[1] $lstym[1] $species[$k]dd1/$unit1a/$unit1c >> pgtmp0
#       endif
# only use Ferret's scaling:
        echo "let vmax=0." >> pgtmp0
       foreach kk($m_pl)
           echo "let vmax1=1.05*$species[$k][de1,$pll1$pld[$kk]$lim2]/$unit1a/$unit1b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
       end
#       echo "if "'`'"vmax eq 0."'`'" then let vmax=1 endif" >> pgtmp0
#       plot
       if ($pl_ty != 3) then
           echo $plcmd$pld[1] $lstym[1] $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
       else
           echo $plcmd$pld[1] $lstym[1] $species[$k]dd1/$unit1a/$unit1c >> pgtmp0
       endif
    endif

# 2 runs
    if ($n_runs == 2) then
#       find max to correctly scale axes
        echo "let vmax=0." >> pgtmp0
       foreach kk($m_pl)
           echo "let vmax1=1.05*$species[$k][de1,$pll1$pld[$kk]$lim2]/$unit1a/$unit1b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de3,$pll1$pld[$kk]$lim2]/$unit3a/$unit3b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
        end
       echo "if "'`'"vmax eq 0."'`'" then let vmax=1 endif" >> pgtmp0
#       plot
       echo $plcmd$pld[1] $lstym[1] $lim1'`'"1.05*vmax"'`' $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
       echo $plcmdo$pld[1] $lstym[1] $xax$species[$k]dd3/$unit3a/$unit3c$yax >> pgtmp0
    endif
# 3 runs
    if ($n_runs == 3) then
#       find max to correctly scale axes
        echo "let vmax=0." >> pgtmp0
       foreach kk($m_pl)
           echo "let vmax1=1.05*$species[$k][de1,$pll1$pld[$kk]$lim2]/$unit1a/$unit1b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de3,$pll1$pld[$kk]$lim2]/$unit3a/$unit3b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de5,$pll1$pld[$kk]$lim2]/$unit5a/$unit5b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
        end
       echo "if "'`'"vmax eq 0."'`'" then let vmax=1 endif" >> pgtmp0
#       plot
       echo $plcmd$pld[1] $lstym[1] $lim1'`'"1.05*vmax"'`' $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
       echo $plcmdo$pld[1] $lstym[1] $xax$species[$k]dd3/$unit3a/$unit3c$yax >> pgtmp0
       echo $plcmdo$pld[1] $lstym[1] $xax$species[$k]dd5/$unit5a/$unit5c$yax >> pgtmp0
    endif
# four runs
    if ($n_runs == 4) then
#       find max to correctly scale axes
        echo "let vmax=0." >> pgtmp0
       foreach kk($m_pl)
           echo "let vmax1=1.05*$species[$k][de1,$pll1$pld[$kk]$lim2]/$unit1a/$unit1b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de3,$pll1$pld[$kk]$lim2]/$unit3a/$unit3b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de5,$pll1$pld[$kk]$lim2]/$unit5a/$unit5b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
           echo "let vmax1=1.05*$species[$k][de7,$pll1$pld[$kk]$lim2]/$unit7a/$unit7b" >> pgtmp0
           echo "if "'`'"vmax1 gt vmax"'`'" then let vmax="'`'"vmax1"'`'" endif" >> pgtmp0
        end
       echo "if "'`'"vmax eq 0."'`'" then let vmax=1 endif" >> pgtmp0
#       plot
       echo $plcmd$pld[1] $lstym[1] $lim1'`'"1.05*vmax"'`' $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
       echo $plcmdo$pld[1] $lstym[1] $xax$species[$k]dd3/$unit3a/$unit3c$yax >> pgtmp0
       echo $plcmdo$pld[1] $lstym[1] $xax$species[$k]dd5/$unit5a/$unit5c$yax >> pgtmp0
       echo $plcmdo$pld[1] $lstym[1] $xax$species[$k]dd7/$unit7a/$unit7c$yax >> pgtmp0
    endif
#   add name of species
    echo "let label_y = 1.15*vmax" >> pgtmp0
    echo "if "'`'"$pl_ty eq 2 "'`'" then let label_y="'`'"1.05*eta[de2,l=1,k=$maxheight]"'`'" endif" >> pgtmp0
#    echo "if "'`'"$pl_ty eq 3 "'`'" then let label_y="'`'"1.05*eta[de2,l=1,k=$maxheight]"'`'" endif" >> pgtmp0
    echo "if "'`'"$pl_ty eq 3 "'`'" then let label_y=1.05*$maxheight endif" >> pgtmp0
    echo "label 0.," '`'"label_y"'`'",0,0,.2 @P1$species[$k]" >> pgtmp0
#   heights/times to be overplotted
    foreach kk ($m_pl)
        if ($n_runs == 1 && $kk != 1) echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
       if ($n_runs == 2 && $kk != 1) then
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd3/$unit3a/$unit3c$yax >> pgtmp0
       endif
       if ($n_runs == 3 && $kk != 1) then
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd3/$unit3a/$unit3c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd5/$unit5a/$unit5c$yax >> pgtmp0
       endif
       if ($n_runs == 4 && $kk != 1) then
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd1/$unit1a/$unit1c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd3/$unit3a/$unit3c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd5/$unit5a/$unit5c$yax >> pgtmp0
           echo $plcmdo$pld[$kk] $lstym[$kk] $xax$species[$k]dd7/$unit7a/$unit7c$yax >> pgtmp0
       endif
    end

    if ($countVP == $pl_np) then
#       next page
       echo "CANCEL MODE METAFILE" >> pgtmp0
       @ countPG+=1
       echo "SET MODE METAFILE:"$pl_name.$countPG".plt" >> pgtmp0
        set countVP = 1
        echo "SET WINDOW /clear" >> pgtmp0
#       put basic info again on VP1
        echo "set viewport VP$countVP" >> pgtmp0
        echo "plot /i=0:100 /hlimits=0:100 /vlimits=0:100 /noaxis /nolabel 0" >> pgtmp0
        echo "label 10,110,  -1, 0, .25 @P1Gas phase" >> pgtmp0
       echo "label 0, 90, -1, 0, .15 @P1$line1" >> pgtmp0
       echo "plot /overlay /vs /nolabel /line=1 $lsty[1] {80,99,99}, {95,95,95}" >> pgtmp0
       echo "label 0, 80, -1, 0, .15 @P1$line2" >> pgtmp0
       if ($countl >= 2) echo "plot /overlay /vs /nolabel /line=2 $lsty[2] {80,99,99}, {85,85,85}" >> pgtmp0
       echo "label 0, 70, -1, 0, .15 @P1$line3" >> pgtmp0
       if ($countl >= 3) echo "plot /overlay /vs /nolabel /line=3 $lsty[3] {80,99,99}, {75,75,75}" >> pgtmp0
       echo "label 0, 60, -1, 0, .15 @P1$line4" >> pgtmp0
       if ($countl >= 4) echo "plot /overlay /vs /nolabel /line=4 $lsty[4] {80,99,99}, {65,65,65}" >> pgtmp0
       echo "label 0, 50, -1, 0, .15 @P1$line5" >> pgtmp0
       if ($countl >= 5) echo "plot /overlay /vs /nolabel /line=5 $lsty[5] {80,99,99}, {55,55,55}" >> pgtmp0
       echo "label 0, 40, -1, 0, .15 @P1$line6" >> pgtmp0
       if ($countl >= 6) echo "plot /overlay /vs /nolabel /line=6 $lsty[6] {80,99,99}, {45,45,45}" >> pgtmp0
       echo "label 0, 30, -1, 0, .15 @P1$line7" >> pgtmp0
       if ($countl >= 7) echo "plot /overlay /vs /nolabel /line=7 $lsty[7] {80,99,99}, {35,35,35}" >> pgtmp0
       echo "label 0, 20, -1, 0, .15 @P1$line8" >> pgtmp0
       if ($countl >= 8) echo "plot /overlay /vs /nolabel /line=8 $lsty[8] {80,99,99}, {25,25,25}" >> pgtmp0
       echo "label 0, 10, -1, 0, .15 @P1$line9" >> pgtmp0
       if ($countl >= 9) echo "plot /overlay /vs /nolabel /line=9 $lsty[9] {80,99,99}, {15,15,15}" >> pgtmp0
    endif

end

echo "CANCEL MODE METAFILE" >> pgtmp0
# if total number of plots is multiple of "plots per page" there is an empty output file: ignore it!
if ($countVP == 1) @ countPG-=1

# finalize plot file
# d1, d3, d5, .. have to be replaced with sed in final ferret file to give [d=1] etc
sed 's/dd1/[d=1]/g' pgtmp0 > pgtmp1
sed 's/dd2/[d=2]/g' pgtmp1 > pgtmp2
sed 's/dd3/[d=3]/g' pgtmp2 > pgtmp3
sed 's/dd4/[d=4]/g' pgtmp3 > pgtmp4
sed 's/dd5/[d=5]/g' pgtmp4 > pgtmp5
sed 's/dd6/[d=6]/g' pgtmp5 > pgtmp6
sed 's/dd7/[d=7]/g' pgtmp6 > pgtmp6a
sed 's/dd8/[d=8]/g' pgtmp6a > pgtmp6b
sed 's/de1/d=1/g' pgtmp6b > pgtmp7
sed 's/de2/d=2/g' pgtmp7 > pgtmp8
sed 's/de3/d=3/g' pgtmp8 > pgtmp9
sed 's/de4/d=4/g' pgtmp9 > pgtmp10
sed 's/de5/d=5/g' pgtmp10 > pgtmp11
sed 's/de6/d=6/g' pgtmp11 > pgtmp12
sed 's/de7/d=7/g' pgtmp12 > pgtmp13
sed 's/de8/d=8/g' pgtmp13 > pgtmp14

# name of plot: pl_name
mv -f pgtmp14 $pl_name.jnl

# plot

# source /soft/ferret_paths_RH9
# # ferret -batch $pl_name.ps -script $pl_name.jnl
ferret -script $pl_name.jnl

# determine number of meta print files
if ($countPG == 1) set metafiles="$pl_name.1.plt"
if ($countPG == 2) set metafiles="$pl_name.1.plt $pl_name.2.plt"
if ($countPG == 3) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt"
if ($countPG == 4) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt"
if ($countPG == 5) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt"
if ($countPG == 6) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt"
if ($countPG == 7) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt"
if ($countPG == 8) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt"
if ($countPG == 9) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt"
if ($countPG == 10) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt  $pl_name.10.plt"
if ($countPG == 11) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt  $pl_name.10.plt $pl_name.11.plt"
if ($countPG == 12) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt  $pl_name.10.plt $pl_name.11.plt $pl_name.12.plt"
if ($countPG == 13) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt  $pl_name.10.plt $pl_name.11.plt $pl_name.12.plt $pl_name.13.plt"
if ($countPG == 58) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt $pl_name.10.plt $pl_name.11.plt $pl_name.12.plt $pl_name.13.plt $pl_name.14.plt $pl_name.15.plt $pl_name.16.plt $pl_name.17.plt $pl_name.18.plt $pl_name.19.plt $pl_name.20.plt $pl_name.21.plt $pl_name.22.plt $pl_name.23.plt $pl_name.24.plt $pl_name.25.plt $pl_name.26.plt $pl_name.27.plt $pl_name.28.plt $pl_name.29.plt $pl_name.30.plt $pl_name.31.plt $pl_name.32.plt $pl_name.33.plt $pl_name.34.plt $pl_name.35.plt $pl_name.36.plt $pl_name.37.plt $pl_name.38.plt $pl_name.39.plt $pl_name.40.plt $pl_name.41.plt $pl_name.42.plt $pl_name.43.plt $pl_name.44.plt $pl_name.45.plt $pl_name.46.plt $pl_name.47.plt $pl_name.48.plt $pl_name.49.plt $pl_name.50.plt $pl_name.51.plt $pl_name.52.plt $pl_name.53.plt $pl_name.54.plt $pl_name.55.plt $pl_name.56.plt $pl_name.57.plt $pl_name.58.plt"
if ($countPG == 59) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt $pl_name.10.plt $pl_name.11.plt $pl_name.12.plt $pl_name.13.plt $pl_name.14.plt $pl_name.15.plt $pl_name.16.plt $pl_name.17.plt $pl_name.18.plt $pl_name.19.plt $pl_name.20.plt $pl_name.21.plt $pl_name.22.plt $pl_name.23.plt $pl_name.24.plt $pl_name.25.plt $pl_name.26.plt $pl_name.27.plt $pl_name.28.plt $pl_name.29.plt $pl_name.30.plt $pl_name.31.plt $pl_name.32.plt $pl_name.33.plt $pl_name.34.plt $pl_name.35.plt $pl_name.36.plt $pl_name.37.plt $pl_name.38.plt $pl_name.39.plt $pl_name.40.plt $pl_name.41.plt $pl_name.42.plt $pl_name.43.plt $pl_name.44.plt $pl_name.45.plt $pl_name.46.plt $pl_name.47.plt $pl_name.48.plt $pl_name.49.plt $pl_name.50.plt $pl_name.51.plt $pl_name.52.plt $pl_name.53.plt $pl_name.54.plt $pl_name.55.plt $pl_name.56.plt $pl_name.57.plt $pl_name.58.plt $pl_name.59.plt"
if ($countPG == 60) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt $pl_name.10.plt $pl_name.11.plt $pl_name.12.plt $pl_name.13.plt $pl_name.14.plt $pl_name.15.plt $pl_name.16.plt $pl_name.17.plt $pl_name.18.plt $pl_name.19.plt $pl_name.20.plt $pl_name.21.plt $pl_name.22.plt $pl_name.23.plt $pl_name.24.plt $pl_name.25.plt $pl_name.26.plt $pl_name.27.plt $pl_name.28.plt $pl_name.29.plt $pl_name.30.plt $pl_name.31.plt $pl_name.32.plt $pl_name.33.plt $pl_name.34.plt $pl_name.35.plt $pl_name.36.plt $pl_name.37.plt $pl_name.38.plt $pl_name.39.plt $pl_name.40.plt $pl_name.41.plt $pl_name.42.plt $pl_name.43.plt $pl_name.44.plt $pl_name.45.plt $pl_name.46.plt $pl_name.47.plt $pl_name.48.plt $pl_name.49.plt $pl_name.50.plt $pl_name.51.plt $pl_name.52.plt $pl_name.53.plt $pl_name.54.plt $pl_name.55.plt $pl_name.56.plt $pl_name.57.plt $pl_name.58.plt $pl_name.59.plt $pl_name.60.plt"
if ($countPG == 61) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt $pl_name.10.plt $pl_name.11.plt $pl_name.12.plt $pl_name.13.plt $pl_name.14.plt $pl_name.15.plt $pl_name.16.plt $pl_name.17.plt $pl_name.18.plt $pl_name.19.plt $pl_name.20.plt $pl_name.21.plt $pl_name.22.plt $pl_name.23.plt $pl_name.24.plt $pl_name.25.plt $pl_name.26.plt $pl_name.27.plt $pl_name.28.plt $pl_name.29.plt $pl_name.30.plt $pl_name.31.plt $pl_name.32.plt $pl_name.33.plt $pl_name.34.plt $pl_name.35.plt $pl_name.36.plt $pl_name.37.plt $pl_name.38.plt $pl_name.39.plt $pl_name.40.plt $pl_name.41.plt $pl_name.42.plt $pl_name.43.plt $pl_name.44.plt $pl_name.45.plt $pl_name.46.plt $pl_name.47.plt $pl_name.48.plt $pl_name.49.plt $pl_name.50.plt $pl_name.51.plt $pl_name.52.plt $pl_name.53.plt $pl_name.54.plt $pl_name.55.plt $pl_name.56.plt $pl_name.57.plt $pl_name.58.plt $pl_name.59.plt $pl_name.60.plt $pl_name.61.plt"
if ($countPG == 63) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt $pl_name.10.plt $pl_name.11.plt $pl_name.12.plt $pl_name.13.plt $pl_name.14.plt $pl_name.15.plt $pl_name.16.plt $pl_name.17.plt $pl_name.18.plt $pl_name.19.plt $pl_name.20.plt $pl_name.21.plt $pl_name.22.plt $pl_name.23.plt $pl_name.24.plt $pl_name.25.plt $pl_name.26.plt $pl_name.27.plt $pl_name.28.plt $pl_name.29.plt $pl_name.30.plt $pl_name.31.plt $pl_name.32.plt $pl_name.33.plt $pl_name.34.plt $pl_name.35.plt $pl_name.36.plt $pl_name.37.plt $pl_name.38.plt $pl_name.39.plt $pl_name.40.plt $pl_name.41.plt $pl_name.42.plt $pl_name.43.plt $pl_name.44.plt $pl_name.45.plt $pl_name.46.plt $pl_name.47.plt $pl_name.48.plt $pl_name.49.plt $pl_name.50.plt $pl_name.51.plt $pl_name.52.plt $pl_name.53.plt $pl_name.54.plt $pl_name.55.plt $pl_name.56.plt $pl_name.57.plt $pl_name.58.plt $pl_name.59.plt $pl_name.60.plt $pl_name.61.plt $pl_name.62.plt $pl_name.63.plt"

# if ($countPG == 100) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt $pl_name.10.plt $pl_name.11.plt $pl_name.12.plt $pl_name.13.plt $pl_name.14.plt $pl_name.15.plt $pl_name.16.plt $pl_name.17.plt $pl_name.18.plt $pl_name.19.plt $pl_name.20.plt $pl_name.21.plt $pl_name.22.plt $pl_name.23.plt $pl_name.24.plt $pl_name.25.plt $pl_name.26.plt $pl_name.27.plt $pl_name.28.plt $pl_name.29.plt $pl_name.30.plt $pl_name.31.plt $pl_name.32.plt $pl_name.33.plt $pl_name.34.plt $pl_name.35.plt $pl_name.36.plt $pl_name.37.plt $pl_name.38.plt $pl_name.39.plt $pl_name.40.plt $pl_name.41.plt $pl_name.42.plt $pl_name.43.plt $pl_name.44.plt $pl_name.45.plt $pl_name.46.plt $pl_name.47.plt $pl_name.48.plt $pl_name.49.plt $pl_name.50.plt $pl_name.51.plt $pl_name.52.plt $pl_name.53.plt $pl_name.54.plt $pl_name.55.plt $pl_name.56.plt $pl_name.57.plt $pl_name.58.plt $pl_name.59.plt $pl_name.60.plt $pl_name.61.plt $pl_name.62.plt $pl_name.63.plt $pl_name.64.plt $pl_name.65.plt $pl_name.66.plt $pl_name.67.plt $pl_name.68.plt $pl_name.69.plt $pl_name.70.plt $pl_name.71.plt $pl_name.72.plt $pl_name.73.plt $pl_name.74.plt $pl_name.75.plt $pl_name.76.plt $pl_name.77.plt $pl_name.78.plt $pl_name.79.plt $pl_name.80.plt $pl_name.81.plt $pl_name.82.plt $pl_name.83.plt $pl_name.84.plt $pl_name.85.plt $pl_name.86.plt $pl_name.87.plt $pl_name.88.plt $pl_name.89.plt $pl_name.90.plt $pl_name.91.plt $pl_name.92.plt $pl_name.93.plt $pl_name.94.plt $pl_name.95.plt $pl_name.96.plt $pl_name.97.plt $pl_name.98.plt $pl_name.99.plt $pl_name.100.plt"
# if ($countPG == 113) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt $pl_name.10.plt $pl_name.11.plt $pl_name.12.plt $pl_name.13.plt $pl_name.14.plt $pl_name.15.plt $pl_name.16.plt $pl_name.17.plt $pl_name.18.plt $pl_name.19.plt $pl_name.20.plt $pl_name.21.plt $pl_name.22.plt $pl_name.23.plt $pl_name.24.plt $pl_name.25.plt $pl_name.26.plt $pl_name.27.plt $pl_name.28.plt $pl_name.29.plt $pl_name.30.plt $pl_name.31.plt $pl_name.32.plt $pl_name.33.plt $pl_name.34.plt $pl_name.35.plt $pl_name.36.plt $pl_name.37.plt $pl_name.38.plt $pl_name.39.plt $pl_name.40.plt $pl_name.41.plt $pl_name.42.plt $pl_name.43.plt $pl_name.44.plt $pl_name.45.plt $pl_name.46.plt $pl_name.47.plt $pl_name.48.plt $pl_name.49.plt $pl_name.50.plt $pl_name.51.plt $pl_name.52.plt $pl_name.53.plt $pl_name.54.plt $pl_name.55.plt $pl_name.56.plt $pl_name.57.plt $pl_name.58.plt $pl_name.59.plt $pl_name.60.plt $pl_name.61.plt $pl_name.62.plt $pl_name.63.plt $pl_name.64.plt $pl_name.65.plt $pl_name.66.plt $pl_name.67.plt $pl_name.68.plt $pl_name.69.plt $pl_name.70.plt $pl_name.71.plt $pl_name.72.plt $pl_name.73.plt $pl_name.74.plt $pl_name.75.plt $pl_name.76.plt $pl_name.77.plt $pl_name.78.plt $pl_name.79.plt $pl_name.80.plt $pl_name.81.plt $pl_name.82.plt $pl_name.83.plt $pl_name.84.plt $pl_name.85.plt $pl_name.86.plt $pl_name.87.plt $pl_name.88.plt $pl_name.89.plt $pl_name.90.plt $pl_name.91.plt $pl_name.92.plt $pl_name.93.plt $pl_name.94.plt $pl_name.95.plt $pl_name.96.plt $pl_name.97.plt $pl_name.98.plt $pl_name.99.plt $pl_name.100.plt $pl_name.101.plt $pl_name.102.plt $pl_name.103.plt $pl_name.104.plt $pl_name.105.plt $pl_name.106.plt $pl_name.107.plt $pl_name.108.plt $pl_name.109.plt $pl_name.110.plt $pl_name.111.plt $pl_name.112.plt $pl_name.113.plt"
# if ($countPG == 114) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt $pl_name.10.plt $pl_name.11.plt $pl_name.12.plt $pl_name.13.plt $pl_name.14.plt $pl_name.15.plt $pl_name.16.plt $pl_name.17.plt $pl_name.18.plt $pl_name.19.plt $pl_name.20.plt $pl_name.21.plt $pl_name.22.plt $pl_name.23.plt $pl_name.24.plt $pl_name.25.plt $pl_name.26.plt $pl_name.27.plt $pl_name.28.plt $pl_name.29.plt $pl_name.30.plt $pl_name.31.plt $pl_name.32.plt $pl_name.33.plt $pl_name.34.plt $pl_name.35.plt $pl_name.36.plt $pl_name.37.plt $pl_name.38.plt $pl_name.39.plt $pl_name.40.plt $pl_name.41.plt $pl_name.42.plt $pl_name.43.plt $pl_name.44.plt $pl_name.45.plt $pl_name.46.plt $pl_name.47.plt $pl_name.48.plt $pl_name.49.plt $pl_name.50.plt $pl_name.51.plt $pl_name.52.plt $pl_name.53.plt $pl_name.54.plt $pl_name.55.plt $pl_name.56.plt $pl_name.57.plt $pl_name.58.plt $pl_name.59.plt $pl_name.60.plt $pl_name.61.plt $pl_name.62.plt $pl_name.63.plt $pl_name.64.plt $pl_name.65.plt $pl_name.66.plt $pl_name.67.plt $pl_name.68.plt $pl_name.69.plt $pl_name.70.plt $pl_name.71.plt $pl_name.72.plt $pl_name.73.plt $pl_name.74.plt $pl_name.75.plt $pl_name.76.plt $pl_name.77.plt $pl_name.78.plt $pl_name.79.plt $pl_name.80.plt $pl_name.81.plt $pl_name.82.plt $pl_name.83.plt $pl_name.84.plt $pl_name.85.plt $pl_name.86.plt $pl_name.87.plt $pl_name.88.plt $pl_name.89.plt $pl_name.90.plt $pl_name.91.plt $pl_name.92.plt $pl_name.93.plt $pl_name.94.plt $pl_name.95.plt $pl_name.96.plt $pl_name.97.plt $pl_name.98.plt $pl_name.99.plt $pl_name.100.plt $pl_name.101.plt $pl_name.102.plt $pl_name.103.plt $pl_name.104.plt $pl_name.105.plt $pl_name.106.plt $pl_name.107.plt $pl_name.108.plt $pl_name.109.plt $pl_name.110.plt $pl_name.111.plt $pl_name.112.plt $pl_name.113.plt $pl_name.114.plt"
# if ($countPG == 115) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt $pl_name.10.plt $pl_name.11.plt $pl_name.12.plt $pl_name.13.plt $pl_name.14.plt $pl_name.15.plt $pl_name.16.plt $pl_name.17.plt $pl_name.18.plt $pl_name.19.plt $pl_name.20.plt $pl_name.21.plt $pl_name.22.plt $pl_name.23.plt $pl_name.24.plt $pl_name.25.plt $pl_name.26.plt $pl_name.27.plt $pl_name.28.plt $pl_name.29.plt $pl_name.30.plt $pl_name.31.plt $pl_name.32.plt $pl_name.33.plt $pl_name.34.plt $pl_name.35.plt $pl_name.36.plt $pl_name.37.plt $pl_name.38.plt $pl_name.39.plt $pl_name.40.plt $pl_name.41.plt $pl_name.42.plt $pl_name.43.plt $pl_name.44.plt $pl_name.45.plt $pl_name.46.plt $pl_name.47.plt $pl_name.48.plt $pl_name.49.plt $pl_name.50.plt $pl_name.51.plt $pl_name.52.plt $pl_name.53.plt $pl_name.54.plt $pl_name.55.plt $pl_name.56.plt $pl_name.57.plt $pl_name.58.plt $pl_name.59.plt $pl_name.60.plt $pl_name.61.plt $pl_name.62.plt $pl_name.63.plt $pl_name.64.plt $pl_name.65.plt $pl_name.66.plt $pl_name.67.plt $pl_name.68.plt $pl_name.69.plt $pl_name.70.plt $pl_name.71.plt $pl_name.72.plt $pl_name.73.plt $pl_name.74.plt $pl_name.75.plt $pl_name.76.plt $pl_name.77.plt $pl_name.78.plt $pl_name.79.plt $pl_name.80.plt $pl_name.81.plt $pl_name.82.plt $pl_name.83.plt $pl_name.84.plt $pl_name.85.plt $pl_name.86.plt $pl_name.87.plt $pl_name.88.plt $pl_name.89.plt $pl_name.90.plt $pl_name.91.plt $pl_name.92.plt $pl_name.93.plt $pl_name.94.plt $pl_name.95.plt $pl_name.96.plt $pl_name.97.plt $pl_name.98.plt $pl_name.99.plt $pl_name.100.plt $pl_name.101.plt $pl_name.102.plt $pl_name.103.plt $pl_name.104.plt $pl_name.105.plt $pl_name.106.plt $pl_name.107.plt $pl_name.108.plt $pl_name.109.plt $pl_name.110.plt $pl_name.111.plt $pl_name.112.plt $pl_name.113.plt $pl_name.114.plt $pl_name.115.plt"
# if ($countPG == 201) set metafiles="$pl_name.1.plt $pl_name.2.plt $pl_name.3.plt $pl_name.4.plt $pl_name.5.plt $pl_name.6.plt $pl_name.7.plt $pl_name.8.plt $pl_name.9.plt $pl_name.10.plt $pl_name.11.plt $pl_name.12.plt $pl_name.13.plt $pl_name.14.plt $pl_name.15.plt $pl_name.16.plt $pl_name.17.plt $pl_name.18.plt $pl_name.19.plt $pl_name.20.plt $pl_name.21.plt $pl_name.22.plt $pl_name.23.plt $pl_name.24.plt $pl_name.25.plt $pl_name.26.plt $pl_name.27.plt $pl_name.28.plt $pl_name.29.plt $pl_name.30.plt $pl_name.31.plt $pl_name.32.plt $pl_name.33.plt $pl_name.34.plt $pl_name.35.plt $pl_name.36.plt $pl_name.37.plt $pl_name.38.plt $pl_name.39.plt $pl_name.40.plt $pl_name.41.plt $pl_name.42.plt $pl_name.43.plt $pl_name.44.plt $pl_name.45.plt $pl_name.46.plt $pl_name.47.plt $pl_name.48.plt $pl_name.49.plt $pl_name.50.plt $pl_name.51.plt $pl_name.52.plt $pl_name.53.plt $pl_name.54.plt $pl_name.55.plt $pl_name.56.plt $pl_name.57.plt $pl_name.58.plt $pl_name.59.plt $pl_name.60.plt $pl_name.61.plt $pl_name.62.plt $pl_name.63.plt $pl_name.64.plt $pl_name.65.plt $pl_name.66.plt $pl_name.67.plt $pl_name.68.plt $pl_name.69.plt $pl_name.70.plt $pl_name.71.plt $pl_name.72.plt $pl_name.73.plt $pl_name.74.plt $pl_name.75.plt $pl_name.76.plt $pl_name.77.plt $pl_name.78.plt $pl_name.79.plt $pl_name.80.plt $pl_name.81.plt $pl_name.82.plt $pl_name.83.plt $pl_name.84.plt $pl_name.85.plt $pl_name.86.plt $pl_name.87.plt $pl_name.88.plt $pl_name.89.plt $pl_name.90.plt $pl_name.91.plt $pl_name.92.plt $pl_name.93.plt $pl_name.94.plt $pl_name.95.plt $pl_name.96.plt $pl_name.97.plt $pl_name.98.plt $pl_name.99.plt $pl_name.100.plt $pl_name.101.plt $pl_name.102.plt $pl_name.103.plt $pl_name.104.plt $pl_name.105.plt $pl_name.106.plt $pl_name.107.plt $pl_name.108.plt $pl_name.109.plt $pl_name.110.plt $pl_name.111.plt $pl_name.112.plt $pl_name.113.plt $pl_name.114.plt $pl_name.115.plt $pl_name.116.plt $pl_name.117.plt $pl_name.118.plt $pl_name.119.plt $pl_name.120.plt $pl_name.121.plt $pl_name.122.plt $pl_name.123.plt $pl_name.124.plt $pl_name.125.plt $pl_name.126.plt $pl_name.127.plt $pl_name.128.plt $pl_name.129.plt $pl_name.130.plt $pl_name.131.plt $pl_name.132.plt $pl_name.133.plt $pl_name.134.plt $pl_name.135.plt $pl_name.136.plt $pl_name.137.plt $pl_name.138.plt $pl_name.139.plt $pl_name.140.plt $pl_name.141.plt $pl_name.142.plt $pl_name.143.plt $pl_name.144.plt $pl_name.145.plt $pl_name.146.plt $pl_name.147.plt $pl_name.148.plt $pl_name.149.plt $pl_name.150.plt $pl_name.151.plt $pl_name.152.plt $pl_name.153.plt $pl_name.154.plt $pl_name.155.plt $pl_name.156.plt $pl_name.157.plt $pl_name.158.plt $pl_name.159.plt $pl_name.160.plt $pl_name.161.plt $pl_name.162.plt $pl_name.163.plt $pl_name.164.plt $pl_name.165.plt $pl_name.166.plt $pl_name.167.plt $pl_name.168.plt $pl_name.169.plt $pl_name.170.plt $pl_name.171.plt $pl_name.172.plt $pl_name.173.plt $pl_name.174.plt $pl_name.175.plt $pl_name.176.plt $pl_name.177.plt $pl_name.178.plt $pl_name.179.plt $pl_name.180.plt $pl_name.181.plt $pl_name.182.plt $pl_name.183.plt $pl_name.184.plt $pl_name.185.plt $pl_name.186.plt $pl_name.187.plt $pl_name.188.plt $pl_name.189.plt $pl_name.190.plt $pl_name.191.plt $pl_name.192.plt $pl_name.193.plt $pl_name.194.plt $pl_name.195.plt $pl_name.196.plt $pl_name.197.plt $pl_name.198.plt $pl_name.199.plt $pl_name.200.plt $pl_name.201.plt"


gksm2ps -p portrait -l cps -d cps -o $pl_name.pre.ps $metafiles

# add page numbering to ps file and delete empty pages
ps2ps $pl_name.pre.ps $pl_name.ps



# clean up ----
# temporary files
rm -rf pgtmp*
# prelim PS file
rm -f $pl_name.pre.ps
# meta print files
rm -f *.plt
