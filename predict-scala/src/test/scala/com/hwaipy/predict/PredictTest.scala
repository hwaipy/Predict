package com.hwaipy.predict

import java.io.File
import org.scalatest._

class PredictTest extends FunSuite with BeforeAndAfter with BeforeAndAfterAll {

  override def beforeAll() {
  }

  override def afterAll() {
  }

  before {
  }

  after {
  }

  test("Test Output.") {
    val expectString = "14051.236329999999725 0.008 148.770 2520.964532688913550\n14051.236341574074686 0.070 148.723 2514.064181249298599\n14051.236353148149647 0.132 148.675 2507.165396287208296\n14051.236364722224607 0.195 148.627 2500.268199907343387\n14051.236376296299568 0.257 148.579 2493.372614429590612\n14051.236387870374529 0.320 148.531 2486.478662391193666\n14051.236399444449489 0.383 148.482 2479.586366550176535\n14051.236411018524450 0.446 148.433 2472.695749887794591\n14051.236422592599411 0.510 148.384 2465.806558482383480\n14051.236434166674371 0.573 148.335 2458.919370100345077\n14051.236445740749332 0.637 148.285 2452.033931212527023\n14051.236457314824293 0.700 148.235 2445.150265724457313\n14051.236468888899253 0.764 148.185 2438.268397780032956\n14051.236480462974214 0.828 148.134 2431.388351764947402\n14051.236492037049175 0.893 148.083 2424.510152310348531\n14051.236503611124135 0.957 148.032 2417.633547675714453\n14051.236515185199096 1.022 147.980 2410.759116308511693\n14051.236526759274057 1.087 147.929 2403.886606899567141\n14051.236538333349017 1.152 147.876 2397.016045093594130\n14051.236549907423978 1.217 147.824 2390.147456797493305\n14051.236561481498939 1.282 147.771 2383.280868184021074\n14051.236573055573899 1.348 147.718 2376.416305694943730\n14051.236584629648860 1.414 147.665 2369.553519984874129\n14051.236596203723821 1.479 147.611 2362.693090249407305\n14051.236607777798781 1.546 147.557 2355.834767617243870\n14051.236619351873742 1.612 147.503 2348.978579644625370\n14051.236630925948702 1.678 147.448 2342.124554175316007\n14051.236642500023663 1.745 147.393 2335.272719345789028\n14051.236654074098624 1.812 147.338 2328.423103588568210\n14051.236665648173584 1.879 147.282 2321.575735635980891\n14051.236677222248545 1.946 147.226 2314.730369170525591\n14051.236688796323506 2.014 147.170 2307.887584340752255\n14051.236700370398466 2.081 147.113 2301.047135357964180\n14051.236711944473427 2.149 147.056 2294.209052197873461\n14051.236723518548388 2.217 146.998 2287.373365157961416\n14051.236735092623348 2.285 146.941 2280.540104861824148\n14051.236746666698309 2.354 146.882 2273.709302263363497\n14051.236758240773270 2.423 146.824 2266.880713976877814\n14051.236769814848230 2.491 146.765 2260.054921085083151\n14051.236781388923191 2.560 146.706 2253.231680786097968\n14051.236792962998152 2.630 146.646 2246.411025403981512\n14051.236804537073112 2.699 146.586 2239.592987617799736\n14051.236816111148073 2.769 146.526 2232.777600466230524\n14051.236827685223034 2.839 146.465 2225.964897352370826\n14051.236839259297994 2.909 146.404 2219.154638116323440\n14051.236850833372955 2.979 146.342 2212.347404884515072\n14051.236862407447916 3.050 146.281 2205.542958143216765\n14051.236873981522876 3.121 146.218 2198.741332807267554\n14051.236885555597837 3.192 146.155 2191.942564183587820\n14051.236897129672798 3.263 146.092 2185.146687976397061\n14051.236908703747758 3.335 146.029 2178.353740292345265\n14051.236920277822719 3.406 145.965 2171.563484523671377\n14051.236931851897680 3.478 145.900 2164.776503969429996\n14051.236943425972640 3.551 145.836 2157.992562736950276\n14051.236955000047601 3.623 145.770 2151.211698603961850\n14051.236966574122562 3.696 145.705 2144.433949781124738\n14051.236978148197522 3.769 145.639 2137.659354918160716\n14051.236989722272483 3.842 145.572 2130.887953110680883\n14051.237001296347444 3.915 145.505 2124.119783906166049\n14051.237012870422404 3.989 145.438 2117.354615203126286\n14051.237024444497365 4.063 145.370 2110.593031823711044\n14051.237036018572326 4.137 145.301 2103.834802471545572\n14051.237047592647286 4.211 145.233 2097.079968570347773\n14051.237059166722247 4.286 145.163 2090.328572029609859\n14051.237070740797208 4.361 145.094 2083.580655253180339\n14051.237082314872168 4.436 145.023 2076.836261145003391\n14051.237093888947129 4.511 144.953 2070.095161983896560\n14051.237105463022090 4.587 144.882 2063.357944111910456\n14051.237117037097050 4.663 144.810 2056.624380702669896\n14051.237128611172011 4.739 144.738 2049.894516739208029\n14051.237140185246972 4.815 144.665 2043.168397743952710\n14051.237151759321932 4.892 144.592 2036.446069786598400\n14051.237163333396893 4.969 144.518 2029.727579491999450\n14051.237174907471854 5.047 144.444 2023.012703976879038\n14051.237186481546814 5.124 144.369 2016.302031306181334\n14051.237198055621775 5.202 144.294 2009.595339592609662\n14051.237209629696736 5.280 144.218 2002.892677767853002\n14051.237221203771696 5.358 144.142 1996.194095362450980\n14051.237232777846657 5.437 144.065 1989.499642514449079\n14051.237244351921618 5.516 143.988 1982.809369978905124\n14051.237255925996578 5.595 143.910 1976.123329136713437\n14051.237267500071539 5.675 143.831 1969.441303262676684\n14051.237279074146500 5.755 143.752 1962.763882679597828\n14051.237290648221460 5.835 143.672 1956.090851789490671\n14051.237302222296421 5.916 143.592 1949.422264575058307\n14051.237313796371382 5.996 143.511 1942.758175695586942\n14051.237325370446342 6.077 143.430 1936.098640494907386\n14051.237336944521303 6.159 143.348 1929.443715013982228\n14051.237348518596264 6.240 143.265 1922.793188533431703\n14051.237360092671224 6.322 143.182 1916.147653645656646\n14051.237371666746185 6.405 143.098 1909.506900889171447\n14051.237383240821146 6.487 143.014 1902.870989199410587\n14051.237394814896106 6.570 142.929 1896.239978263940657\n14051.237406388971067 6.654 142.843 1889.613928533292892\n14051.237417963046028 6.737 142.756 1882.992901233406883\n14051.237429537120988 6.821 142.669 1876.376692299534625\n14051.237441111195949 6.906 142.582 1869.765896910069841\n14051.237452685270910 6.990 142.493 1863.160312405343120\n14051.237464259345870 7.075 142.404 1856.560003235761769\n14051.237475833420831 7.160 142.315 1849.965034688767901\n14051.237487407495792 7.246 142.224 1843.375472902448564\n14051.237498981570752 7.332 142.133 1836.791384878964209\n14051.237510555645713 7.418 142.041 1830.212838498661540\n14051.237522129720674 7.505 141.949 1823.639638194982581\n14051.237533703795634 7.592 141.856 1817.072382549451731\n14051.237545277870595 7.679 141.762 1810.510877594299245\n14051.237556851945556 7.767 141.667 1803.955194859679295\n14051.237568426020516 7.855 141.572 1797.405406825111868\n14051.237580000095477 7.944 141.475 1790.861586933994204\n14051.237591574170438 8.033 141.379 1784.323809608552438\n14051.237603148245398 8.122 141.281 1777.791887599341635\n14051.237614722320359 8.211 141.182 1771.266422917326736\n14051.237626296395320 8.301 141.083 1764.747230102283083\n14051.237637870470280 8.392 140.983 1758.234387650526060\n14051.237649444545241 8.482 140.882 1751.727975117542655\n14051.237661018620202 8.574 140.781 1745.228073135112709\n14051.237672592695162 8.665 140.678 1738.734763428102951\n14051.237684166770123 8.757 140.575 1732.247867988244934\n14051.237695740845083 8.849 140.471 1725.767992740096133\n14051.237707314920044 8.942 140.366 1719.294961680204096\n14051.237718888995005 9.035 140.260 1712.828861082663252\n14051.237730463069965 9.129 140.153 1706.369778405665784\n14051.237742037144926 9.223 140.045 1699.917802308912314\n14051.237753611219887 9.317 139.937 1693.473022673274272\n14051.237765185294847 9.412 139.828 1687.035530620856889\n14051.237776759369808 9.507 139.717 1680.605159979188556\n14051.237788333444769 9.602 139.606 1674.182521824656305\n14051.237799907519729 9.698 139.494 1667.767452267888302\n14051.237811481594690 9.795 139.381 1661.360047596033610\n14051.237823055669651 9.892 139.267 1654.960405439503575\n14051.237834629744611 9.989 139.152 1648.568624793539584\n14051.237846203819572 10.087 139.036 1642.184806040176454\n14051.237857777894533 10.185 138.919 1635.808794616337082\n14051.237869351969493 10.283 138.801 1629.441206783646066\n14051.237880926044454 10.383 138.682 1623.081890539047890\n14051.237892500119415 10.482 138.562 1616.730952037249608\n14051.237904074194375 10.582 138.441 1610.388498935117241\n14051.237915648269336 10.683 138.318 1604.054640414608002\n14051.237927222344297 10.783 138.195 1597.729487208050841\n14051.237938796419257 10.885 138.071 1591.412897674032820\n14051.237950370494218 10.987 137.946 1585.105493978088816\n14051.237961944569179 11.089 137.819 1578.807137347503158\n14051.237973518644139 11.192 137.692 1572.517944966615914\n14051.237985092719100 11.295 137.563 1566.238035697526357\n14051.237996666794061 11.399 137.433 1559.967530107738639\n14051.238008240869021 11.503 137.302 1553.706550496001682\n14051.238019814943982 11.608 137.170 1547.455220920703596\n14051.238031389018943 11.713 137.036 1541.213416308255319\n14051.238042963093903 11.819 136.902 1534.981766562015991\n14051.238054537168864 11.925 136.766 1528.760149873173532\n14051.238066111243825 12.032 136.629 1522.548697622558166\n14051.238077685318785 12.139 136.491 1516.347543095451556\n14051.238089259393746 12.246 136.351 1510.156821510622194\n14051.238100833468707 12.355 136.211 1503.976670051301198\n14051.238112407543667 12.463 136.068 1497.806979896202620\n14051.238123981618628 12.573 135.925 1491.648388690193087\n14051.238135555693589 12.682 135.780 1485.500791263771589\n14051.238147129768549 12.793 135.634 1479.364332972138300\n14051.238158703843510 12.904 135.487 1473.239161293892494\n14051.238170277918471 13.015 135.338 1467.125425864122235\n14051.238181851993431 13.127 135.188 1461.023278508953808\n14051.238193426068392 13.239 135.036 1454.932628478903325\n14051.238205000143353 13.352 134.883 1448.854122162690146\n14051.238216574218313 13.466 134.729 1442.787672887839108\n14051.238228148293274 13.580 134.573 1436.733441585497303\n14051.238239722368235 13.695 134.415 1430.691591551174952\n14051.238251296443195 13.810 134.256 1424.662288481210908\n14051.238262870518156 13.926 134.096 1418.645700508985328\n14051.238274444593117 14.042 133.934 1412.641756954775019\n14051.238286018668077 14.159 133.770 1406.651114041571418\n14051.238297592743038 14.276 133.605 1400.673705628787957\n14051.238309166817999 14.394 133.438 1394.709709972401697\n14051.238320740892959 14.513 133.270 1388.759307954760743\n14051.238332314967920 14.632 133.100 1382.822683123383513\n14051.238343889042881 14.751 132.928 1376.900021729686387\n14051.238355463117841 14.872 132.755 1370.991512769646079\n14051.238367037192802 14.993 132.579 1365.097111171967526\n14051.238378611267763 15.114 132.403 1359.217485832141620\n14051.238390185342723 15.236 132.224 1353.352596787812217\n14051.238401759417684 15.359 132.044 1347.502644415078521\n14051.238413333492645 15.482 131.861 1341.667832037661356\n14051.238424907567605 15.605 131.677 1335.848365968662165\n14051.238436481642566 15.730 131.491 1330.044455552322461\n14051.238448055717527 15.855 131.304 1324.256080649877958\n14051.238459629792487 15.980 131.114 1318.483922554716855\n14051.238471203867448 16.106 130.922 1312.727966763791073\n14051.238482777942409 16.233 130.729 1306.988435170071398\n14051.238494352017369 16.360 130.533 1301.265552910859014\n14051.238505926092330 16.488 130.336 1295.559548410679781\n14051.238517500167291 16.617 130.136 1289.870653424228749\n14051.238529074242251 16.746 129.935 1284.198875247168644\n14051.238540648317212 16.876 129.731 1278.544908804293300\n14051.238552222392173 17.006 129.525 1272.908767566307233\n14051.238563796467133 17.137 129.317 1267.290697043681348\n14051.238575370542094 17.268 129.107 1261.690946294314472\n14051.238586944617055 17.400 128.895 1256.109767966283471\n14051.238598518692015 17.533 128.681 1250.547418341576531\n14051.238610092766976 17.666 128.464 1245.004157378513582\n14051.238621666841937 17.800 128.245 1239.480026902314648\n14051.238633240916897 17.934 128.024 1233.975738852063159\n14051.238644814991858 18.069 127.800 1228.491341834578179\n14051.238656389066819 18.205 127.574 1223.027110938793385\n14051.238667963141779 18.341 127.345 1217.583325140661600\n14051.238679537216740 18.478 127.115 1212.160267343895384\n14051.238691111291701 18.615 126.881 1206.758224420181705\n14051.238702685366661 18.753 126.645 1201.377271196803576\n14051.238714259441622 18.891 126.407 1196.018135580578928\n14051.238725833516582 19.030 126.166 1190.680899667474705\n14051.238737407591543 19.170 125.923 1185.365866592007251\n14051.238748981666504 19.310 125.676 1180.073343648092987\n14051.238760555741464 19.451 125.428 1174.803642324087377\n14051.238772129816425 19.592 125.176 1169.557078337955545\n14051.238783703891386 19.734 124.922 1164.333762003257334\n14051.238795277966346 19.876 124.665 1159.134437892472761\n14051.238806852041307 20.019 124.405 1153.959223989633074\n14051.238818426116268 20.162 124.143 1148.808453251197761\n14051.238830000191228 20.306 123.877 1143.682463016233214\n14051.238841574266189 20.450 123.609 1138.581595032883342\n14051.238853148341150 20.595 123.337 1133.506195482888415\n14051.238864722416110 20.740 123.063 1128.456615004202604\n14051.238876296491071 20.886 122.786 1123.433007134277659\n14051.238887870566032 21.032 122.505 1118.436135711876432\n14051.238899444640992 21.179 122.222 1113.466162221358672\n14051.238911018715953 21.326 121.935 1108.523455309373730\n14051.238922592790914 21.474 121.646 1103.608388159701690\n14051.238934166865874 21.621 121.353 1098.721338501765331\n14051.238945740940835 21.770 121.056 1093.862688616456808\n14051.238957315015796 21.918 120.757 1089.032631602303809\n14051.238968889090756 22.068 120.454 1084.231947505968719\n14051.238980463165717 22.217 120.148 1079.460837366357737\n14051.238992037240678 22.367 119.839 1074.719701670506538\n14051.239003611315638 22.517 119.526 1070.008945440069283\n14051.239015185390599 22.667 119.209 1065.328978218366728\n14051.239026759465560 22.818 118.889 1060.680214051100165\n14051.239038333540520 22.969 118.566 1056.062886344724802\n14051.239049907615481 23.120 118.239 1051.477789622843829\n14051.239061481690442 23.271 117.908 1046.925164891946679\n14051.239073055765402 23.423 117.574 1042.405443950449353\n14051.239084629840363 23.574 117.236 1037.919062947745715\n14051.239096203915324 23.726 116.894 1033.466462338544488\n14051.239107777990284 23.878 116.549 1029.048086833694015\n14051.239119352065245 24.031 116.199 1024.664385344962056\n14051.239130926140206 24.183 115.846 1020.315636679360637\n14051.239142500215166 24.335 115.489 1016.002647893591302\n14051.239154074290127 24.487 115.128 1011.725704447866974\n14051.239165648365088 24.639 114.763 1007.485271393535299\n14051.239177222440048 24.792 114.394 1003.281817620496213\n14051.239188796515009 24.944 114.021 999.115815766509400\n14051.239200370589970 25.096 113.644 994.987742119849941\n14051.239211944664930 25.248 113.263 990.897912752517072\n14051.239223518739891 25.399 112.878 986.847140030726109\n14051.239235092814852 25.551 112.488 982.835745221894172\n14051.239246666889812 25.702 112.095 978.864218100074481\n14051.239258240964773 25.853 111.697 974.933051495177438\n14051.239269815039734 26.004 111.295 971.042741152858071\n14051.239281389114694 26.154 110.889 967.193785585002274\n14051.239292963189655 26.304 110.478 963.386533589270812\n14051.239304537264616 26.453 110.063 959.621795093059063\n14051.239316111339576 26.602 109.644 955.899921915987079\n14051.239327685414537 26.750 109.220 952.221421985277971\n14051.239339259489498 26.898 108.792 948.586805127930347\n14051.239350833564458 27.045 108.359 944.996582872926410\n14051.239362407639419 27.192 107.922 941.451268243359777\n14051.239373981714380 27.338 107.481 937.951375541647280\n14051.239385555789340 27.483 107.035 934.497282091526017\n14051.239397129864301 27.627 106.584 931.089782017411494\n14051.239408703939262 27.770 106.130 927.729252185817700\n14051.239420278014222 27.912 105.670 924.416209684732735\n14051.239431852089183 28.054 105.206 921.151171672471605\n14051.239443426164144 28.194 104.738 917.934655109797063\n14051.239455000239104 28.333 104.265 914.767176483723233\n14051.239466574314065 28.471 103.787 911.649127082370342\n14051.239478148389026 28.608 103.305 908.581272489245407\n14051.239489722463986 28.744 102.819 905.563999583927398\n14051.239501296538947 28.878 102.328 902.597820029179275\n14051.239512870613908 29.011 101.832 899.683243508912710\n14051.239524444688868 29.142 101.332 896.820777403461534\n14051.239536018763829 29.272 100.827 894.010926458112408\n14051.239547592838790 29.401 100.318 891.254082606625957\n14051.239559166913750 29.527 99.805 888.550966142167908\n14051.239570740988711 29.652 99.287 885.901959860041529\n14051.239582315063672 29.775 98.765 883.307554510470482\n14051.239593889138632 29.897 98.239 880.768236450755921\n14051.239605463213593 30.016 97.708 878.284487276608047\n14051.239617037288554 30.134 97.173 875.856783449263503\n14051.239628611363514 30.249 96.633 873.485595917626370\n14051.239640185438475 30.363 96.090 871.171297781534463\n14051.239651759513436 30.474 95.542 868.914534051293799\n14051.239663333588396 30.583 94.990 866.715662583436256\n14051.239674907663357 30.690 94.435 864.575128439609898\n14051.239686481738318 30.794 93.875 862.493369239551953\n14051.239698055813278 30.896 93.311 860.470814770073048\n14051.239709629888239 30.996 92.744 858.507886595279047\n14051.239721203963200 31.093 92.173 856.604922320478181\n14051.239732778038160 31.187 91.598 854.762479029838232\n14051.239744352113121 31.279 91.019 852.980873518357043\n14051.239755926188082 31.368 90.437 851.260490582615262\n14051.239767500263042 31.454 89.852 849.601704872013670\n14051.239779074338003 31.537 89.263 848.004880514206206\n14051.239790648412963 31.618 88.671 846.470370746286676\n14051.239802222487924 31.695 88.076 844.998459600028468\n14051.239813796562885 31.770 87.478 843.589595896693481\n14051.239825370637845 31.841 86.877 842.244037574465324\n14051.239836944712806 31.910 86.273 840.962090756145471\n14051.239848518787767 31.975 85.667 839.744048938733613\n14051.239860092862727 32.037 85.058 838.590192674092577\n14051.239871666937688 32.096 84.446 837.500789260070292\n14051.239883241012649 32.151 83.832 836.476052522008331\n14051.239894815087609 32.203 83.216 835.516304833271874\n14051.239906389162570 32.252 82.598 834.621729468312537\n14051.239917963237531 32.297 81.978 833.792537868121826\n14051.239929537312491 32.339 81.357 833.028926866826737\n14051.239941111387452 32.378 80.733 832.331078464212283\n14051.239952685462413 32.412 80.108 831.699159614202699\n14051.239964259537373 32.444 79.482 831.133322029032684\n14051.239975833612334 32.472 78.854 830.633683233277793\n14051.239987407687295 32.496 78.226 830.200404144590266\n14051.239998981762255 32.516 77.596 829.833568314628110\n14051.240010555837216 32.533 76.966 829.533264879836338\n14051.240022129912177 32.547 76.335 829.299567020536870\n14051.240033703987137 32.556 75.704 829.132531873571793\n14051.240045278062098 32.562 75.072 829.032200463314666\n14051.240056852137059 32.565 74.440 828.998597643322682\n14051.240068426212019 32.564 73.808 829.031734789254074\n14051.240080000286980 32.559 73.176 829.131601676257219\n14051.240091574361941 32.550 72.545 829.298174581633020\n14051.240103148436901 32.538 71.914 829.531413619537034\n14051.240114722511862 32.522 71.284 829.831262788781828\n14051.240126296586823 32.503 70.654 830.197650040727808\n14051.240137870661783 32.480 70.025 830.630506115817070\n14051.240149444736744 32.453 69.398 831.129692316574165\n14051.240161018811705 32.423 68.771 831.695105124970269\n14051.240172592886665 32.389 68.146 832.326609355868186\n14051.240184166961626 32.352 67.522 833.024054395974645\n14051.240195741036587 32.311 66.900 833.787274381044995\n14051.240207315111547 32.267 66.280 834.616088390647064\n14051.240218889186508 32.219 65.662 835.510300659506697\n14051.240230463261469 32.168 65.046 836.469740711973486\n14051.240242037336429 32.114 64.432 837.494106580207472\n14051.240253611411390 32.056 63.820 838.583196674007809\n14051.240265185486351 31.995 63.211 839.736758263912179\n14051.240276759561311 31.931 62.604 840.954525055095019\n14051.240288333636272 31.864 62.000 842.236217483992164\n14051.240299907711233 31.793 61.398 843.581543025926635\n14051.240311481786193 31.720 60.800 844.990254459020889\n14051.240323055861154 31.643 60.205 846.461920943547511\n14051.240334629936115 31.564 59.612 847.996268429340375\n14051.240346204011075 31.482 59.023 849.592955836604233\n14051.240357778086036 31.396 58.438 851.251630814403711\n14051.240369352160997 31.308 57.855 852.971930103091950\n14051.240380926235957 31.218 57.277 854.753479902420509\n14051.240392500310918 31.124 56.701 856.595971591303737\n14051.240404074385879 31.028 56.130 858.498863148694113\n14051.240415648460839 30.930 55.562 860.461824318676008\n14051.240427222535800 30.829 54.998 862.484442920644369\n14051.240438796610761 30.725 54.438 864.566298141697757\n14051.240450370685721 30.619 53.882 866.706960926079319\n14051.240461944760682 30.511 53.330 868.905994365179367\n14051.240473518835643 30.401 52.781 871.162954087776257\n14051.240485092910603 30.288 52.238 873.477482917956081\n14051.240496666985564 30.174 51.698 875.848936477122038\n14051.240508241060525 30.057 51.162 878.276942303869873\n14051.240519815135485 29.938 50.631 880.761030059328277\n14051.240531389210446 29.818 50.104 883.300723872243225\n14051.240542963285407 29.695 49.582 885.895542716664977\n14051.240554537360367 29.571 49.063 888.545000785559864\n14051.240566111435328 29.445 48.549 891.248717717311479\n14051.240577685510289 29.318 48.040 894.005981676094734\n14051.240589259585249 29.189 47.535 896.816402391889596\n14051.240600833660210 29.058 47.035 899.679478566691273\n14051.240612407735171 28.926 46.539 902.594705909580512\n14051.240623981810131 28.792 46.047 905.561577476004231\n14051.240635555885092 28.658 45.560 908.579584000069417\n14051.240647129960053 28.521 45.078 911.648338691757544\n14051.240658704035013 28.384 44.600 914.767081669132153\n14051.240670278109974 28.246 44.127 917.935421066534673\n14051.240681852184935 28.106 43.658 921.152841503769992\n14051.240693426259895 27.965 43.194 924.418826825385509\n14051.240705000334856 27.824 42.734 927.732860386353423\n14051.240716574409817 27.681 42.279 931.094425328395118\n14051.240728148484777 27.537 41.828 934.503004848036767\n14051.240739722559738 27.393 41.382 937.958222394114273\n14051.240751296634699 27.248 40.940 941.459284006714142\n14051.240762870709659 27.102 40.503 945.005812646944378\n14051.240774444784620 26.955 40.070 948.597294239006146\n14051.240786018859581 26.808 39.641 952.233215972628727\n14051.240797592934541 26.660 39.217 955.913066519131462\n14051.240809167009502 26.512 38.798 959.636336238400986\n14051.240820741084462 26.363 38.383 963.402669760810795\n14051.240832315159423 26.214 37.972 967.211258342783026\n14051.240843889234384 26.064 37.565 971.061749235995762\n14051.240855463309344 25.914 37.163 974.953641459397204\n14051.240867037384305 25.764 36.765 978.886436630334742\n14051.240878611459266 25.613 36.371 982.859639122127874\n14051.240890185534226 25.462 35.981 986.872756213221578\n14051.240901759609187 25.311 35.596 990.925462062890460\n14051.240913333684148 25.159 35.214 995.016944065548159\n14051.240924907759108 25.008 34.837 999.146881291003297\n14051.240936481834069 24.856 34.464 1003.314793977495810\n14051.240948055909030 24.704 34.095 1007.520205900643987\n14051.240959629983990 24.552 33.730 1011.762644477958133\n14051.240971204058951 24.401 33.369 1016.041640867805086\n14051.240982778133912 24.249 33.012 1020.356730059535607\n14051.240994352208872 24.097 32.658 1024.707626713639911\n14051.241005926283833 23.945 32.309 1029.093523627557943\n14051.241017500358794 23.794 31.963 1033.514142100394338\n14051.241029074433754 23.642 31.621 1037.969033234457811\n14051.241040648508715 23.491 31.283 1042.457752327631624\n14051.241052222583676 23.340 30.949 1046.979858927687474\n14051.241063796658636 23.189 30.618 1051.534916883251753\n14051.241075370733597 23.038 30.291 1056.122679610326941\n14051.241086944808558 22.887 29.967 1060.742350526531709\n14051.241098518883518 22.737 29.647 1065.393690554222076\n14051.241110092958479 22.587 29.331 1070.076281144798031\n14051.241121667033440 22.437 29.018 1074.789708227016035\n14051.241133241108400 22.288 28.709 1079.533562228735946\n14051.241144815183361 22.139 28.402 1084.307438095381713\n14051.241156389258322 21.990 28.100 1089.111129155225171\n14051.241167963333282 21.842 27.800 1093.943852891933602\n14051.241179537408243 21.694 27.504 1098.805410540255480\n14051.241191111483204 21.546 27.211 1103.695415219316828\n14051.241202685558164 21.399 26.921 1108.613484600032052\n14051.241214259633125 21.252 26.635 1113.559240902738566\n14051.241225833708086 21.106 26.352 1118.532310890878989\n14051.241237407783046 20.960 26.071 1123.532325862588550\n14051.241248981858007 20.814 25.794 1128.559124406078809\n14051.241260555932968 20.669 25.520 1133.611942366880385\n14051.241272130007928 20.525 25.249 1138.690626274509441\n14051.241283704082889 20.381 24.980 1143.794825428621834\n14051.241295278157850 20.237 24.715 1148.924193584003433\n14051.241306852232810 20.094 24.453 1154.078388927604692\n14051.241318426307771 19.951 24.193 1159.257074053873112\n14051.241330000382732 19.809 23.936 1164.460125748701785\n14051.241341574457692 19.668 23.682 1169.686796676117183\n14051.241353148532653 19.527 23.431 1174.936971335472663\n14051.241364722607614 19.386 23.182 1180.210329676419178\n14051.241376296682574 19.246 22.936 1185.506555910437783\n14051.241387870757535 19.107 22.693 1190.825338477784953\n14051.241399444832496 18.968 22.452 1196.166370013024334\n14051.241411018907456 18.830 22.214 1201.529563517384076\n14051.241422592982417 18.692 21.979 1206.914188355523265\n14051.241434167057378 18.555 21.746 1212.320164858956332\n14051.241445741132338 18.419 21.515 1217.747202075002861\n14051.241457315207299 18.283 21.287 1223.195013058790892\n14051.241468889282260 18.147 21.061 1228.663314833337836\n14051.241480463357220 18.012 20.838 1234.151828349346260\n14051.241492037432181 17.878 20.617 1239.660500465690575\n14051.241503611507142 17.745 20.398 1245.188616609921382\n14051.241515185582102 17.612 20.181 1250.736130493840619\n14051.241526759657063 17.479 19.967 1256.302778367275778\n14051.241538333732024 17.347 19.755 1261.888300199060495\n14051.241549907806984 17.216 19.546 1267.492439635110486\n14051.241561481881945 17.085 19.338 1273.114943955395574\n14051.241573055956906 16.955 19.133 1278.755564030307141\n14051.241584630031866 16.825 18.929 1284.414282294384293\n14051.241596204106827 16.697 18.728 1290.090401342194582\n14051.241607778181788 16.568 18.529 1295.783909859846290\n14051.241619352256748 16.440 18.332 1301.494572648586654\n14051.241630926331709 16.313 18.137 1307.222157884148146\n14051.241642500406670 16.187 17.944 1312.966437073337602\n14051.241654074481630 16.061 17.753 1318.727185010726998\n14051.241665648556591 15.935 17.563 1324.504412487389800\n14051.241677222631552 15.811 17.376 1330.297435881440833\n14051.241688796706512 15.686 17.191 1336.106271695717851\n14051.241700370781473 15.563 17.007 1341.930707446777660\n14051.241711944856434 15.440 16.825 1347.770533725192308\n14051.241723518931394 15.317 16.645 1353.625544153876717\n14051.241735093006355 15.195 16.467 1359.495535345105509\n14051.241746667081316 15.074 16.291 1365.380543918522108\n14051.241758241156276 14.953 16.116 1371.279898807262271\n14051.241769815231237 14.833 15.943 1377.193641813388467\n14051.241781389306198 14.713 15.772 1383.121581096241698\n14051.241792963381158 14.594 15.603 1389.063527599640338\n14051.241804537456119 14.476 15.435 1395.019295012901921\n14051.241816111531080 14.358 15.268 1400.988699729570726\n14051.241827685606040 14.241 15.104 1406.971560809649191\n14051.241839259681001 14.124 14.941 1412.967941449385535\n14051.241850833755962 14.008 14.779 1418.977183428223498\n14051.241862407830922 13.892 14.619 1424.999354549905092\n14051.241873981905883 13.777 14.461 1431.034284152310192\n14051.241885555980843 13.663 14.304 1437.081804047484866\n14051.241897130055804 13.549 14.149 1443.141748485921426\n14051.241908704130765 13.435 13.995 1449.213954119380787\n14051.241920278205725 13.322 13.842 1455.298504998243743\n14051.241931852280686 13.210 13.691 1461.394752880036776\n14051.241943426355647 13.098 13.542 1467.502785955169429\n14051.241955000430607 12.987 13.394 1473.622450124445550\n14051.241966574505568 12.876 13.247 1479.753593513230044\n14051.241978148580529 12.766 13.101 1485.896066438040407\n14051.241989722655489 12.656 12.957 1492.049721372999329\n14051.242001296730450 12.547 12.815 1498.214661162651055\n14051.242012870805411 12.438 12.673 1504.390246439705152\n14051.242024444880371 12.330 12.533 1510.576583758773040\n14051.242036018955332 12.222 12.394 1516.773533897008065\n14051.242047593030293 12.115 12.257 1522.980959627108632\n14051.242059167105253 12.008 12.120 1529.198725687489969\n14051.242070741180214 11.902 11.985 1535.426698750783089\n14051.242082315255175 11.797 11.852 1541.664747394546794\n14051.242093889330135 11.692 11.719 1547.912993646795940\n14051.242105463405096 11.587 11.587 1554.170807045627271\n14051.242117037480057 11.483 11.457 1560.438312883855360\n14051.242128611555017 11.379 11.328 1566.715387064215747\n14051.242140185629978 11.276 11.200 1573.001907249565647\n14051.242151759704939 11.173 11.073 1579.297752835439496\n14051.242163333779899 11.071 10.948 1585.602804922803671\n14051.242174907854860 10.969 10.823 1591.917200511791634\n14051.242186481929821 10.868 10.700 1598.240315952321680\n14051.242198056004781 10.767 10.577 1604.572291159377755\n14051.242209630079742 10.666 10.456 1610.913013791876438\n14051.242221204154703 10.567 10.336 1617.262373084470482\n14051.242232778229663 10.467 10.216 1623.620259821893114\n14051.242244352304624 10.368 10.098 1629.986566315218852\n14051.242255926379585 10.269 9.981 1636.361443016226303\n14051.242267500454545 10.171 9.865 1642.744272268382929\n14051.242279074529506 10.074 9.750 1649.135207124283397\n14051.242290648604467 9.976 9.635 1655.534145763167089\n14051.242302222679427 9.879 9.522 1661.940987773990855\n14051.242313796754388 9.783 9.410 1668.355634132479963\n14051.242325370829349 9.687 9.299 1674.777987180303626\n14051.242336944904309 9.592 9.188 1681.207950602931078\n14051.242348518979270 9.496 9.079 1687.645688559983000\n14051.242360093054231 9.402 8.970 1694.090589355851307\n14051.242371667129191 9.307 8.863 1700.542819432595934\n14051.242383241204152 9.214 8.756 1707.002287643327236\n14051.242394815279113 9.120 8.650 1713.468904082112658\n14051.242406389354073 9.027 8.545 1719.942580064701133\n14051.242417963429034 8.934 8.441 1726.423228109784304\n14051.242429537503995 8.842 8.338 1732.911023072060971\n14051.242441111578955 8.750 8.235 1739.405357785660499\n14051.242452685653916 8.659 8.134 1745.906409141781523\n14051.242464259728877 8.568 8.033 1752.414094283348049\n14051.242475833803837 8.477 7.933 1758.928331463975383\n14051.242487407878798 8.387 7.834 1765.449040031904133\n14051.242498981953759 8.297 7.736 1771.976140410737571\n14051.242510556028719 8.207 7.638 1778.509817071488214\n14051.242522130103680 8.118 7.542 1785.049466816094764\n14051.242533704178641 8.029 7.446 1791.595275933340872\n14051.242545278253601 7.941 7.351 1798.147168987490431\n14051.242556852328562 7.852 7.256 1804.705071538231778\n14051.242568426403523 7.765 7.162 1811.268910124796776\n14051.242580000478483 7.677 7.070 1817.838612250622191\n14051.242591574553444 7.590 6.977 1824.414106368397597\n14051.242603148628405 7.504 6.886 1830.995586764037853\n14051.242614722703365 7.417 6.795 1837.582454173836140\n14051.242626296778326 7.331 6.705 1844.174904480607665\n14051.242637870853287 7.246 6.616 1850.772869789448350\n14051.242649444928247 7.161 6.527 1857.376283083521230\n14051.242661019003208 7.076 6.439 1863.985078210223264\n14051.242672593078169 6.991 6.352 1870.599189867857376\n14051.242684167153129 6.907 6.265 1877.218820017344569\n14051.242695741228090 6.823 6.179 1883.843372379772063\n14051.242707315303051 6.739 6.094 1890.473050343193336\n14051.242718889378011 6.656 6.009 1897.107791881750700\n14051.242730463452972 6.573 5.925 1903.747535757362584\n14051.242742037527933 6.490 5.842 1910.392221507484919\n14051.242753611602893 6.408 5.759 1917.041789433509848\n14051.242765185677854 6.326 5.677 1923.696448413552162\n14051.242776759752815 6.244 5.595 1930.355604782979754\n14051.242788333827775 6.163 5.514 1937.019468697081720\n14051.242799907902736 6.082 5.434 1943.687983398158622\n14051.242811481977697 6.001 5.354 1950.361092835279578\n14051.242823056052657 5.921 5.275 1957.038741654648902\n14051.242834630127618 5.841 5.197 1963.720875188666469\n14051.242846204202579 5.761 5.119 1970.407708556785565\n14051.242857778277539 5.681 5.041 1977.098650385181372\n14051.242869352352500 5.602 4.965 1983.793916937128188\n14051.242880926427461 5.523 4.888 1990.493456188115488\n14051.242892500502421 5.444 4.813 1997.197216749364770\n14051.242904074577382 5.366 4.738 2003.905147858178225\n14051.242915648652342 5.288 4.663 2010.617199369957007\n14051.242927222727303 5.210 4.589 2017.333321747372793\n14051.242938796802264 5.132 4.515 2024.053736506157293\n14051.242950370877224 5.055 4.442 2030.777854548622599\n14051.242961944952185 4.978 4.370 2037.505898400923570\n14051.242973519027146 4.901 4.298 2044.237820868577728\n14051.242985093102106 4.825 4.226 2050.973575320996133\n14051.242996667177067 4.749 4.156 2057.713115683432079\n14051.243008241252028 4.673 4.085 2064.456396428935477\n14051.243019815326988 4.597 4.015 2071.203644096374319\n14051.243031389401949 4.521 3.946 2077.954271321765191\n14051.243042963476910 4.446 3.877 2084.708505551743656\n14051.243054537551870 4.371 3.808 2091.466303369878915\n14051.243066111626831 4.297 3.740 2098.227621868225469\n14051.243077685701792 4.222 3.673 2104.992418640015785\n14051.243089259776752 4.148 3.606 2111.760651772047822\n14051.243100833851713 4.074 3.539 2118.532552350918195\n14051.243112407926674 4.001 3.473 2125.307534537612810\n14051.243123982001634 3.927 3.407 2132.085830234862442\n14051.243135556076595 3.854 3.342 2138.867399433893297\n14051.243147130151556 3.781 3.277 2145.652202584927181\n14051.243158704226516 3.709 3.213 2152.440200590695895\n14051.243170278301477 3.636 3.149 2159.231354799343535\n14051.243181852376438 3.564 3.085 2166.025626998966345\n14051.243193426451398 3.492 3.022 2172.823252952039184\n14051.243205000526359 3.420 2.960 2179.623648344800586\n14051.243216574601320 3.349 2.897 2186.427049664576316\n14051.243228148676280 3.277 2.836 2193.233420394855330\n14051.243239722751241 3.206 2.774 2200.042724428206384\n14051.243251296826202 3.136 2.713 2206.854926060226262\n14051.243262870901162 3.065 2.653 2213.669989983729465\n14051.243274444976123 2.995 2.592 2220.488155646602536\n14051.243286019051084 2.924 2.533 2227.308839905567311\n14051.243297593126044 2.855 2.473 2234.132282862135526\n14051.243309167201005 2.785 2.414 2240.958450743853064\n14051.243320741275966 2.715 2.356 2247.787310147624794\n14051.243332315350926 2.646 2.297 2254.618828036008836\n14051.243343889425887 2.577 2.240 2261.452971731248454\n14051.243355463500848 2.508 2.182 2268.289984026664115\n14051.243367037575808 2.439 2.125 2275.129282818295906\n14051.243378611650769 2.371 2.068 2281.971111491285228\n14051.243390185725730 2.303 2.012 2288.815438757368611\n14051.243401759800690 2.235 1.956 2295.662233664313590\n14051.243413333875651 2.167 1.900 2302.511465590664102\n14051.243424907950612 2.099 1.845 2309.363104240906978\n14051.243436482025572 2.032 1.790 2316.217119641670706\n14051.243448056100533 1.964 1.736 2323.073758038250617\n14051.243459630175494 1.897 1.681 2329.932438377685230\n14051.243471204250454 1.830 1.627 2336.793407433012362\n14051.243482778325415 1.764 1.574 2343.656636474332572\n14051.243494352400376 1.697 1.521 2350.522097072464931\n14051.243505926475336 1.631 1.468 2357.389761094051664\n14051.243517500550297 1.565 1.415 2364.259600698320355\n14051.243529074625258 1.499 1.363 2371.131864858396966\n14051.243540648700218 1.433 1.311 2378.005973338764761\n14051.243552222775179 1.367 1.260 2384.882175592121257\n14051.243563796850140 1.302 1.208 2391.760444907044985\n14051.243575370925100 1.237 1.157 2398.640754845752781\n14051.243586945000061 1.172 1.107 2405.523079239833805\n14051.243598519075022 1.107 1.057 2412.407392186268225\n14051.243610093149982 1.042 1.007 2419.293945141643690\n14051.243621667224943 0.978 0.957 2426.182158607225574\n14051.243633241299904 0.914 0.907 2433.072284474744720\n14051.243644815374864 0.849 0.858 2439.964297869849361\n14051.243656389449825 0.785 0.810 2446.858174165995933\n14051.243667963524786 0.722 0.761 2453.753888982050739\n14051.243679537599746 0.658 0.713 2460.651418179144912\n14051.243691111674707 0.594 0.665 2467.550737857121931\n14051.243702685749668 0.531 0.617 2474.452102040128011\n14051.243714259824628 0.468 0.570 2481.354931987997588\n14051.243725833899589 0.405 0.523 2488.259482116623531\n14051.243737407974550 0.342 0.476 2495.165729449941864\n14051.243748982049510 0.280 0.430 2502.073651235537454\n14051.243760556124471 0.217 0.384 2508.983224941730896\n14051.243772130199432 0.155 0.338 2515.894428254535796\n14051.243783704274392 0.093 0.292 2522.807517232437021\n14051.243795278349353 0.031 0.247 2529.721913736620081\n14051.243806852424314 -0.031 0.201 2536.637874183217264\n14051.243818426499274 -0.093 0.157 2543.555377101708473\n14051.243830000574235 -0.155 0.112 2550.474401224944813\n14051.243841574649196 -0.216 0.068 2557.394925488191348\n14051.243853148724156 -0.277 0.024 2564.316929025141690\n14051.243864722799117 -0.338 359.980 2571.240669748633536\n14051.243876296874078 -0.399 359.936 2578.165570073877461\n14051.243887870949038 -0.460 359.893 2585.091888241050128"
    val expect = new PassDetail
    expectString.split("\n").toList.foreach(line => {
      val items = line.split(" +")
      expect.appendPosition(items(0).toDouble, items(1).toDouble, items(2).toDouble, items(3).toDouble)
    })
    val predict = new Predict
    predict.ReadDataFiles(new File("predict.tle"))
    val actual = predict.Predict(Array[String]("predict", "14051.23633", "32.326", "80.026", "5075.0", "1"))
    assert(actual == expect)
  }
}