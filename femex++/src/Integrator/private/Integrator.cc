/*
 * Integrator.h
 *
 *  Created on: Oct 9, 2014
 *      Author: lurker
 */


#include "Integrator.h"

void Integrator::clear() {

	qpts.clear();
	qwts.clear();

}

void Integrator::GaussData(){

	if ( prec > 23 ){
		mexErrMsgTxt("Precision No More than 23\n");
	}

	Real_t nodes[]={
	   0.0000000000000000e+00,
	  -5.7735026918962573e-01,
	   5.7735026918962573e-01,
	  -7.7459666924148340e-01,
	   0.0000000000000000e+00,
	   7.7459666924148340e-01,
	  -8.6113631159405302e-01,
	  -3.3998104358485626e-01,
	   3.3998104358485626e-01,
	   8.6113631159405291e-01,
	  -9.0617984593866396e-01,
	  -5.3846931010568289e-01,
	   0.0000000000000000e+00,
	   5.3846931010568300e-01,
	   9.0617984593866385e-01,
	  -9.3246951420315205e-01,
	  -6.6120938646626459e-01,
	  -2.3861918608319677e-01,
	   2.3861918608319677e-01,
	   6.6120938646626470e-01,
	   9.3246951420315205e-01,
	  -9.4910791234275904e-01,
	  -7.4153118559939368e-01,
	  -4.0584515137739718e-01,
	   0.0000000000000000e+00,
	   4.0584515137739707e-01,
	   7.4153118559939446e-01,
	   9.4910791234275904e-01,
	  -9.6028985649753695e-01,
	  -7.9666647741362617e-01,
	  -5.2553240991632888e-01,
	  -1.8343464249564978e-01,
	   1.8343464249564989e-01,
	   5.2553240991632899e-01,
	   7.9666647741362606e-01,
	   9.6028985649753684e-01,
	  -9.6816023950762597e-01,
	  -8.3603110732663699e-01,
	  -6.1337143270059025e-01,
	  -3.2425342340380880e-01,
	   0.0000000000000000e+00,
	   3.2425342340380880e-01,
	   6.1337143270059080e-01,
	   8.3603110732663488e-01,
	   9.6816023950762697e-01,
	  -9.7390652851717197e-01,
	  -8.6506336668898398e-01,
	  -6.7940956829902455e-01,
	  -4.3339539412924699e-01,
	  -1.4887433898163116e-01,
	   1.4887433898163116e-01,
	   4.3339539412924699e-01,
	   6.7940956829902444e-01,
	   8.6506336668898398e-01,
	   9.7390652851717197e-01,
	  -9.7822865814603999e-01,
	  -8.8706259976811996e-01,
	  -7.3015200657404225e-01,
	  -5.1909612920681159e-01,
	  -2.6954315595234490e-01,
	   0.0000000000000000e+00,
	   2.6954315595234490e-01,
	   5.1909612920681170e-01,
	   7.3015200657405002e-01,
	   8.8706259976809299e-01,
	   9.7822865814605797e-01,
	  -9.8156063424673201e-01,
	  -9.0411725637045204e-01,
	  -7.6990267419431768e-01,
	  -5.8731795428661426e-01,
	  -3.6783149899818035e-01,
	  -1.2523340851146880e-01,
	   1.2523340851146880e-01,
	   3.6783149899818035e-01,
	   5.8731795428661426e-01,
	   7.6990267419431779e-01,
	   9.0411725637045204e-01,
	   9.8156063424673201e-01};

	Real_t weights[]={
	   2.0000000000000000e+00,
	   1.0000000000000000e+00,
	   1.0000000000000000e+00,
	   5.5555555555555525e-01,
	   8.8888888888888895e-01,
	   5.5555555555555525e-01,
	   3.4785484513745391e-01,
	   6.5214515486254621e-01,
	   6.5214515486254621e-01,
	   3.4785484513745391e-01,
	   2.3692688505618870e-01,
	   4.7862867049936653e-01,
	   5.6888888888888889e-01,
	   4.7862867049936653e-01,
	   2.3692688505618870e-01,
	   1.7132449237917091e-01,
	   3.6076157304813788e-01,
	   4.6791393457269131e-01,
	   4.6791393457269131e-01,
	   3.6076157304813788e-01,
	   1.7132449237917091e-01,
	   1.2948496616886801e-01,
	   2.7970539148927831e-01,
	   3.8183005050511859e-01,
	   4.1795918367346940e-01,
	   3.8183005050511881e-01,
	   2.7970539148927598e-01,
	   1.2948496616886970e-01,
	   1.0122853629037380e-01,
	   2.2238103445337859e-01,
	   3.1370664587788738e-01,
	   3.6268378337836188e-01,
	   3.6268378337836188e-01,
	   3.1370664587788738e-01,
	   2.2238103445337859e-01,
	   1.0122853629037380e-01,
	   8.1274388361575897e-02,
	   1.8064816069485429e-01,
	   2.6061069640293560e-01,
	   3.1234707704000292e-01,
	   3.3023935500125973e-01,
	   3.1234707704000247e-01,
	   2.6061069640293533e-01,
	   1.8064816069485770e-01,
	   8.1274388361572095e-02,
	   6.6671344308686806e-02,
	   1.4945134915057301e-01,
	   2.1908636251598321e-01,
	   2.6926671930999679e-01,
	   2.9552422471475293e-01,
	   2.9552422471475293e-01,
	   2.6926671930999679e-01,
	   2.1908636251598321e-01,
	   1.4945134915057301e-01,
	   6.6671344308686806e-02,
	   5.5668567116215838e-02,
	   1.2558036946487430e-01,
	   1.8629021092774040e-01,
	   2.3319376459199270e-01,
	   2.6280454451024660e-01,
	   2.7292508677790062e-01,
	   2.6280454451024660e-01,
	   2.3319376459199331e-01,
	   1.8629021092773390e-01,
	   1.2558036946491319e-01,
	   5.5668567116169583e-02,
	   4.7175336386475468e-02,
	   1.0693932599536370e-01,
	   1.6007832854335860e-01,
	   2.0316742672306720e-01,
	   2.3349253653835339e-01,
	   2.4914704581340269e-01,
	   2.4914704581340269e-01,
	   2.3349253653835339e-01,
	   2.0316742672306720e-01,
	   1.6007832854335860e-01,
	   1.0693932599536370e-01,
	   4.7175336386475468e-02};

	std::size_t n=floor((prec +2)/2);

	std::size_t ptr1=(n-1)*n/2+1;
	std::size_t ptr2=n*(n+1)/2;
	qpts.resize(ptr2 - ptr1 + 1);
	qwts.resize(ptr2 - ptr1 + 1);

	for (std::size_t _index = 0; _index < qpts.size(); _index++){
		qpts[_index] = nodes[ptr1 - 1 + _index];
		qwts[_index] = weights[ptr1 - 1 + _index];
	}
}

/*
 * All decimal data are with 16 digits significant.
 */
void Integrator::QuadratureData(){

	if (prec > 20) {
		mexErrMsgTxt("Precision No More than 20\n");
	}

	size_t pptr[]={1,2,3,5,7,10,13,17,22,28,34,41,49,59,69,80,93,108,125,142};

	size_t n[]={1,3,1,3,3,3,1,3,3,3,3,6,1,3,3,6,1,3,3,3,6,1,3,3,3,
	3,6,1,3,3,6,6,6,3,3,3,3,3,6,6,3,3,3,3,3,6,6,6,1,3,
	3,3,3,3,3,6,6,6,3,3,3,3,3,3,6,6,6,6,3,3,3,3,3,3,6,
	6,6,6,6,1,3,3,3,3,3,3,3,6,6,6,6,6,1,3,3,3,3,3,3,3,
	3,6,6,6,6,6,6,1,3,3,3,3,3,3,3,3,3,6,6,6,6,6,6,6,1,3,3,3,3,3,3,3,3,
	6,6,6,6,6,6,6,6,1,3,3,3,3,3,3,3,3,3,3,6,6,6,6,6,6,6,6};

	Real_t weights[]={
	1.000000000000000,0.333333333333333,-0.562500000000000,
	0.520833333333333,0.223381589678011,0.109951743655322,
	0.225000000000000,0.132394152788506,0.125939180544827,
	0.116786275726379,0.050844906370207,0.082851075618374,
	-0.149570044467682,0.175615257433208,0.053347235608838,
	0.077113760890257,0.144315607677787,0.095091634267285,
	0.103217370534718,0.032458497623198,0.027230314174435,
	0.097135796282799,0.031334700227139,0.077827541004774,
	0.079647738927210,0.025577675658698,0.043283539377289,
	0.090817990382754,0.036725957756467,0.045321059435528,
	0.072757916845420,0.028327242531057,0.009421666963733,
	0.000927006328961,0.077149534914813,0.059322977380774,
	0.036184540503418,0.013659731002678,0.052337111962204,
	0.020707659639141,0.025731066440455,0.043692544538038,
	0.062858224217885,0.034796112930709,0.006166261051559,
	0.040371557766381,0.022356773202303,0.017316231108659,
	0.052520923400802,0.011280145209330,0.031423518362454,
	0.047072502504194,0.047363586536355,0.031167529045794,
	0.007975771465074,0.036848402728732,0.017401463303822,
	0.015521786839045,0.021883581369429,0.032788353544125,
	0.051774104507292,0.042162588736993,0.014433699669777,
	0.004923403602400,0.024665753212564,0.038571510787061,
	0.014436308113534,0.005010228838501,0.001916875642849,
	0.044249027271145,0.051186548718852,0.023687735870688,
	0.013289775690021,0.004748916608192,0.038550072599593,
	0.027215814320624,0.002182077366797,0.021505319847731,
	0.007673942631049,0.046875697427642,0.006405878578585,
	0.041710296739387,0.026891484250064,0.042132522761650,
	0.030000266842773,0.014200098925024,0.003582462351273,
	0.032773147460627,0.015298306248441,0.002386244192839,
	0.019084792755899,0.006850054546542,0.033437199290803,
	0.005093415440507,0.014670864527638,0.024350878353672,
	0.031107550868969,0.031257111218620,0.024815654339665,
	0.014056073070557,0.003194676173779,0.008119655318993,
	0.026805742283163,0.018459993210822,0.008476868534328,
	0.018292796770025,0.006665632004165,
	0.030809939937647,0.009072436679404,0.018761316939594,
	0.019441097985477,0.027753948610810,0.032256225351457,
	0.025074032616922,0.015271927971832,0.006793922022963,
	-0.002223098729920,0.006331914076406,0.027257538049138,
	0.017676785649465,0.018379484638070,0.008104732808192,
	0.007634129070725,0.000046187660794,0.032906331388919,
	0.010330731891272,0.022387247263016,0.030266125869468,
	0.030490967802198,0.024159212741641,0.016050803586801,
	0.008084580261784,0.002079352027485,0.003884876904981,
	0.025574160612022,0.008880903573338,0.016124546761731,
	0.002491941817491,0.018242840118951,0.010258563736199,
	0.003799928855302,
	0.033057055541624,0.000867019185663,0.011660052716448,
	0.022876936356421,0.030448982673938,0.030624891725355,
	0.024368057676800,0.015997432032024,0.007698301815602,
	-0.000632060497488,0.001751134301193,0.016465839189576,
	0.004839033540485,0.025840906534650,0.008471091054441,
	0.018354914106280,0.000704404677908,0.010112684927462,
	0.003537909385950};
	Real_t coord1[]={
	0.333333333333333,0.666666666666667,0.333333333333333,
	0.600000000000000,0.108103018168070,0.816847572980459,
	0.333333333333333,0.059715871789770,0.797426985353087,
	0.501426509658179,0.873821971016996,0.053145049844817,
	0.333333333333333,0.479308067841920,0.869739794195568,
	0.048690315425316,0.333333333333333,0.081414823414554,
	0.658861384496480,0.898905543365938,0.008394777409958,
	0.333333333333333,0.020634961602525,0.125820817014127,
	0.623592928761935,0.910540973211095,0.036838412054736,
	0.333333333333333,0.028844733232685,0.781036849029926,
	0.141707219414880,0.025003534762686,0.009540815400299,
	-0.069222096541517,0.202061394068290,0.593380199137435,
	0.761298175434837,0.935270103777448,0.050178138310495,
	0.021022016536166,0.023565220452390,0.120551215411079,
	0.457579229975768,0.744847708916828,0.957365299093579,
	0.115343494534698,0.022838332222257,0.025734050548330,
	0.333333333333333,0.009903630120591,0.062566729780852,
	0.170957326397447,0.541200855914337,0.771151009607340,
	0.950377217273082,0.094853828379579,0.018100773278807,
	0.022233076674090,0.022072179275643,0.164710561319092,
	0.453044943382323,0.645588935174913,0.876400233818255,
	0.961218077502598,0.057124757403648,0.092916249356972,
	0.014646950055654,0.001268330932872,-0.013945833716486,
	0.137187291433955,0.444612710305711,0.747070217917492,
	0.858383228050628,0.962069659517853,0.133734161966621,
	0.036366677396917,-0.010174883126571,0.036843869875878,
	0.012459809331199,0.333333333333333,0.005238916103123,
	0.173061122901295,0.059082801866017,0.518892500060958,
	0.704068411554854,0.849069624685052,0.966807194753950,
	0.103575692245252,0.020083411655416,-0.004341002614139,
	0.041941786468010,0.014317320230681,0.333333333333333,
	0.005658918886452,0.035647354750751,0.099520061958437,
	0.199467521245206,0.495717464058095,0.675905990683077,
	0.848248235478508,0.968690546064356,0.010186928826919,
	0.135440871671036,0.054423924290583,0.012868560833637,
	0.067165782413524,0.014663182224828,
	0.333333333333333,0.013310382738157,0.061578811516086,
	0.127437208225989,0.210307658653168,0.500410862393686,
	0.677135612512315,0.846803545029257,0.951495121293100,
	0.913707265566071,0.008430536202420,0.131186551737188,
	0.050203151565675,0.066329263810916,0.011996194566236,
	0.014858100590125,-0.035222015287949,0.333333333333333,
	0.020780025853987,0.090926214604215,0.197166638701138,
	0.488896691193805,0.645844115695741,0.779877893544096,
	0.888942751496321,0.974756272445543,0.003611417848412,
	0.134466754530780,0.014446025776115,0.046933578838178,
	0.002861120350567,0.223861424097916,0.034647074816760,
	0.010161119296270,0.333333333333333,-0.001900928704400,
	0.023574084130543,0.089726636099435,0.196007481363421,
	0.488214180481157,0.647023488009788,0.791658289326483,
	0.893862072318140,0.916762569607942,0.976836157186356,
	0.048741583664839,0.006314115948605,0.134316520547348,
	0.013973893962392,0.075549132909764,-0.008363153208227,
	0.026686063258714,0.010547719294141};
	Real_t coord2[]={
	0.333333333333333,0.166666666666667,0.333333333333333,
	0.200000000000000,0.445948490915965,0.091576213509771,
	0.333333333333333,0.470142064105115,0.101286507323456,
	0.249286745170910,0.063089014491502,0.310352451033784,
	0.333333333333333,0.260345966079040,0.065130102902216,
	0.312865496004874,0.333333333333333,0.459292588292723,
	0.170569307751760,0.050547228317031,0.263112829634638,
	0.333333333333333,0.489682519198738,0.437089591492937,
	0.188203535619033,0.044729513394453,0.221962989160766,
	0.333333333333333,0.485577633383657,0.109481575485037,
	0.307939838764121,0.246672560639903,0.066803251012200,
	0.534611048270758,0.398969302965855,0.203309900431282,
	0.119350912282581,0.032364948111276,0.356620648261293,
	0.171488980304042,0.488217389773805,0.439724392294460,
	0.271210385012116,0.127576145541586,0.021317350453210,
	0.275713269685514,0.281325580989940,0.116251915907597,
	0.333333333333333,0.495048184939705,0.468716635109574,
	0.414521336801277,0.229399572042831,0.114424495196330,
	0.024811391363459,0.268794997058761,0.291730066734288,
	0.126357385491669,0.488963910362179,0.417644719340454,
	0.273477528308839,0.177205532412543,0.061799883090873,
	0.019390961248701,0.172266687821356,0.336861459796345,
	0.298372882136258,0.118974497696957,0.506972916858243,
	0.431406354283023,0.277693644847144,0.126464891041254,
	0.070808385974686,0.018965170241073,0.261311371140087,
	0.388046767090269,0.285712220049916,0.215599664072284,
	0.103575616576386,0.333333333333333,0.497380541948438,
	0.413469438549352,0.470458599066991,0.240553749969521,
	0.147965794222573,0.075465187657474,0.016596402623025,
	0.296555596579887,0.337723063403079,0.204748281642812,
	0.189358492130623,0.085283615682657,0.333333333333333,
	0.497170540556774,0.482176322624625,0.450239969020782,
	0.400266239377397,0.252141267970953,0.162047004658461,
	0.075875882260746,0.015654726967822,0.334319867363658,
	0.292221537796944,0.319574885423190,0.190704224192292,
	0.180483211648746,0.080711313679564,
	0.333333333333333,0.493344808630921,0.469210594241957,
	0.436281395887006,0.394846170673416,0.249794568803157,
	0.161432193743843,0.076598227485371,0.024252439353450,
	0.043146367216965,0.358911494940944,0.294402476751957,
	0.325017801641814,0.184737559666046,0.218796800013321,
	0.101179597136408,0.020874755282586,0.333333333333333,
	0.489609987073006,0.454536892697893,0.401416680649431,
	0.255551654403098,0.177077942152130,0.110061053227952,
	0.055528624251840,0.012621863777229,0.395754787356943,
	0.307929983880436,0.264566948406520,0.358539352205951,
	0.157807405968595,0.075050596975911,0.142421601113383,
	0.065494628082938,0.333333333333333,0.500950464352200,
	0.488212957934729,0.455136681950283,0.401996259318289,
	0.255892909759421,0.176488255995106,0.104170855337758,
	0.053068963840930,0.041618715196029,0.011581921406822,
	0.344855770229001,0.377843269594854,0.306635479062357,
	0.249419362774742,0.212775724802802,0.146965436053239,
	0.137726978828923,0.059696109149007};
	Real_t coord3[]={
	0.333333333333333,0.166666666666667,0.333333333333333,
	0.200000000000000,0.445948490915965,0.091576213509771,
	0.333333333333333,0.470142064105115,0.101286507323456,
	0.249286745170910,0.063089014491502,0.636502499121399,
	0.333333333333333,0.260345966079040,0.065130102902216,
	0.638444188569810,0.333333333333333,0.459292588292723,
	0.170569307751760,0.050547228317031,0.728492392955404,
	0.333333333333333,0.489682519198738,0.437089591492937,
	0.188203535619033,0.044729513394453,0.741198598784498,
	0.333333333333333,0.485577633383657,0.109481575485037,
	0.550352941820999,0.728323904597411,0.923655933587500,
	0.534611048270758,0.398969302965855,0.203309900431282,
	0.119350912282581,0.032364948111276,0.593201213428213,
	0.807489003159792,0.488217389773805,0.439724392294460,
	0.271210385012116,0.127576145541586,0.021317350453210,
	0.608943235779788,0.695836086787803,0.858014033544073,
	0.333333333333333,0.495048184939705,0.468716635109574,
	0.414521336801277,0.229399572042831,0.114424495196330,
	0.024811391363459,0.636351174561660,0.690169159986905,
	0.851409537834241,0.488963910362179,0.417644719340454,
	0.273477528308839,0.177205532412543,0.061799883090873,
	0.019390961248701,0.770608554774996,0.570222290846683,
	0.686980167808088,0.879757171370171,0.506972916858243,
	0.431406354283023,0.277693644847144,0.126464891041254,
	0.070808385974686,0.018965170241073,0.604954466893291,
	0.575586555512814,0.724462663076655,0.747556466051838,
	0.883964574092416,0.333333333333333,0.497380541948438,
	0.413469438549352,0.470458599066991,0.240553749969521,
	0.147965794222573,0.075465187657474,0.016596402623025,
	0.599868711174861,0.642193524941505,0.799592720971327,
	0.768699721401368,0.900399064086661,0.333333333333333,
	0.497170540556774,0.482176322624625,0.450239969020782,
	0.400266239377397,0.252141267970953,0.162047004658461,
	0.075875882260746,0.015654726967822,0.655493203809423,
	0.572337590532020,0.626001190286228,0.796427214974071,
	0.752351005937729,0.904625504095608,
	0.333333333333333,0.493344808630921,0.469210594241957,
	0.436281395887006,0.394846170673416,0.249794568803157,
	0.161432193743843,0.076598227485371,0.024252439353450,
	0.043146367216965,0.632657968856636,0.574410971510855,
	0.624779046792512,0.748933176523037,0.769207005420443,
	0.883962302273467,1.014347260005363,0.333333333333333,
	0.489609987073006,0.454536892697893,0.401416680649431,
	0.255551654403098,0.177077942152130,0.110061053227952,
	0.055528624251840,0.012621863777229,0.600633794794645,
	0.557603261588784,0.720987025817365,0.594527068955871,
	0.839331473680839,0.701087978926173,0.822931324069857,
	0.924344252620784,0.333333333333333,0.500950464352200,
	0.488212957934729,0.455136681950283,0.401996259318289,
	0.255892909759421,0.176488255995106,0.104170855336758,
	0.053068963840930,0.041618715196029,0.011581921406822,
	0.606402646106160,0.615842614456541,0.559048000390295,
	0.736606743262866,0.711675142287434,0.861402717154987,
	0.835586957912363,0.929756171556853};


	size_t i = pptr[prec - 1];
	size_t k, N, m;
	Real_t x, y, z;

	if (prec < 20)
		k = pptr[prec] - i;
	else
		k = sizeof(weights)/sizeof(Real_t) - i + 1;

	N = 0;
	for (size_t j = i - 1; j < i + k - 1 ; j++){
		N += n[j];
	}



	qwts.resize(N);
	vector<vector<Real_t>> bcoords;
	bcoords.resize(N);
	for (size_t j = 0; j < N; j++){
		bcoords[j].resize(3);
	}

	m = 0;


	for (size_t j = i - 1; j < i + k - 1; j++){

		x = coord1[j];
		y = coord2[j];
		z = coord3[j];

		if (n[j] == 1){
			qwts[m] = weights[j];
			bcoords[m][0] = x;
			bcoords[m][1] = y;
			bcoords[m][2] = z;
			m++;
		}
		else if (n[j] == 3){
			qwts[m] = weights[j];
			bcoords[m][0] = x;
			bcoords[m][1] = y;
			bcoords[m][2] = y;
			m++;

			qwts[m] = weights[j];
			bcoords[m][0] = y;
			bcoords[m][1] = x;
			bcoords[m][2] = y;
			m++;

			qwts[m] = weights[j];
			bcoords[m][0] = y;
			bcoords[m][1] = y;
			bcoords[m][2] = x;
			m++;
 		}
		else if (n[j] == 6){
			qwts[m] = weights[j];
			bcoords[m][0] = x;
			bcoords[m][1] = y;
			bcoords[m][2] = z;
			m++;

			qwts[m] = weights[j];
			bcoords[m][0] = x;
			bcoords[m][1] = z;
			bcoords[m][2] = y;
			m++;

			qwts[m] = weights[j];
			bcoords[m][0] = y;
			bcoords[m][1] = x;
			bcoords[m][2] = z;
			m++;

			qwts[m] = weights[j];
			bcoords[m][0] = y;
			bcoords[m][1] = z;
			bcoords[m][2] = x;
			m++;

			qwts[m] = weights[j];
			bcoords[m][0] = z;
			bcoords[m][1] = x;
			bcoords[m][2] = y;
			m++;

			qwts[m] = weights[j];
			bcoords[m][0] = z;
			bcoords[m][1] = y;
			bcoords[m][2] = x;
			m++;
		}
	}

	qpts.resize(2*N);
	Real_t vert[3][2] = { {0, 0}, {1, 0}, {0, 1}};


	for (size_t j = 0; j < N; j++){
		qpts[2*j] = bcoords[j][0]*vert[0][0] +
				bcoords[j][1]*vert[1][0] +
				bcoords[j][2]*vert[2][0];
		qpts[2*j + 1] = bcoords[j][0]*vert[0][1] +
				bcoords[j][1]*vert[1][1] +
				bcoords[j][2]*vert[2][1];

	}

}


namespace {

// Create a new instance of Mesh and return its session id.
MEX_DEFINE(new) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 2);
	OutputArguments output(nlhs, plhs, 1);
	output.set(0, Session<Integrator>::create(new Integrator(input.get<int>(0), input.get<int>(1))));
}

// Delete the instance by id
MEX_DEFINE(delete) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	Session<Integrator>::destroy(input.get(0));
}


MEX_DEFINE(export)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 2);
	Integrator *integrator = Session<Integrator>::get(input.get(0));
	plhs[0] = mxCreateNumericMatrix(integrator->_dim, integrator->qpts.size()/integrator->_dim, mxDOUBLE_CLASS, mxREAL);
	memcpy(mxGetPr(plhs[0]), &(integrator->qpts[0]), integrator->qpts.size()*sizeof(REAL));
	plhs[1] = mxCreateNumericMatrix(integrator->qpts.size()/integrator->_dim, 1,  mxDOUBLE_CLASS, mxREAL);
	memcpy(mxGetPr(plhs[1]), &(integrator->qwts[0]), integrator->qwts.size()*sizeof(REAL));

}

} // namespace


MEX_DISPATCH




