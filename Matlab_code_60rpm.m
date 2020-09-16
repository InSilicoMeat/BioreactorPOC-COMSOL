function out = model
%
% Matlab_code_rpm.m
%
% Model exported on Sep 14 2020, 13:16 by COMSOL 5.5.0.359.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('D:\Jiajia Chen\Cell-based Meat');

model.label('100mL bioreactor_tetra mesh_60rpm_for exporting code.mph');

model.param.set('radius', '5.5[cm]/2');
model.param.set('volume', '100[ml]+volume_impeller');
model.param.set('height', 'volume/(pi*radius^2)');
model.param.set('location_impeller', '1/3*height');
model.param.set('stir_len', '3.8[cm]');
model.param.set('stir_dia', '0.8[cm]');
model.param.set('TIME', '0[s]');
model.param.set('volume_impeller', 'pi*(stir_dia/2)^2*stir_len');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 3);

model.component('comp1').sorder('quadratic');

model.component('comp1').mesh.create('mesh1');
model.component('comp1').mesh.create('mesh2');

model.geom.load({'part1'}, 'Mixer_Module\Impellers,_Axial\pitched_blade_impeller_bent_blade.mph', {'part1'});
model.geom.load({'part2'}, 'Mixer_Module\Impellers,_Radial\backswept_impeller.mph', {'part1'});
model.component('comp1').geom('geom1').create('cyl1', 'Cylinder');
model.component('comp1').geom('geom1').feature('cyl1').set('r', 'radius');
model.component('comp1').geom('geom1').feature('cyl1').set('h', 'height');
model.component('comp1').geom('geom1').create('cyl2', 'Cylinder');
model.component('comp1').geom('geom1').feature('cyl2').set('pos', {'-stir_len/2' '0' 'location_impeller'});
model.component('comp1').geom('geom1').feature('cyl2').set('axis', [1 0 0]);
model.component('comp1').geom('geom1').feature('cyl2').set('r', 'stir_dia/2');
model.component('comp1').geom('geom1').feature('cyl2').set('h', 'stir_len');
model.component('comp1').geom('geom1').create('cyl3', 'Cylinder');
model.component('comp1').geom('geom1').feature('cyl3').set('pos', {'0' '0' 'location_impeller-(stir_dia/2+1.5[cm])/2'});
model.component('comp1').geom('geom1').feature('cyl3').set('axis', [0 0 1]);
model.component('comp1').geom('geom1').feature('cyl3').set('r', '(stir_len)/2+0.5[cm]');
model.component('comp1').geom('geom1').feature('cyl3').set('h', 'stir_dia/2+1.5[cm]');
model.component('comp1').geom('geom1').create('dif1', 'Difference');
model.component('comp1').geom('geom1').feature('dif1').selection('input').set({'cyl1' 'cyl3'});
model.component('comp1').geom('geom1').feature('dif1').selection('input2').set({'cyl2'});
model.component('comp1').geom('geom1').run;

model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material('mat1').propertyGroup('def').func.create('eta', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('Cp', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('rho', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('k', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('cs', 'Interpolation');
model.component('comp1').material('mat1').propertyGroup('def').func.create('an1', 'Analytic');
model.component('comp1').material('mat1').propertyGroup('def').func.create('an2', 'Analytic');
model.component('comp1').material('mat1').propertyGroup('def').func.create('an3', 'Analytic');

model.component('comp1').common.create('rot1', 'RotatingDomain');
model.component('comp1').common('rot1').selection.set([2]);

model.component('comp1').physics.create('spf', 'TurbulentFlowkeps', 'geom1');
model.component('comp1').physics('spf').create('vf1', 'VolumeForce', 3);
model.component('comp1').physics('spf').feature('vf1').selection.set([1 2]);
model.component('comp1').physics('spf').create('sfs1', 'StationaryFreeSurface', 2);
model.component('comp1').physics('spf').feature('sfs1').selection.set([4]);

model.component('comp1').mesh('mesh1').create('ftet1', 'FreeTet');
model.component('comp1').mesh('mesh1').feature('ftet1').create('size1', 'Size');
model.component('comp1').mesh('mesh2').autoMeshSize(2);

model.component('comp1').view('view1').set('renderwireframe', true);

model.component('comp1').material('mat1').label('Water, liquid');
model.component('comp1').material('mat1').set('family', 'water');
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('pieces', {'273.15' '413.15' '1.3799566804-0.021224019151*T^1+1.3604562827E-4*T^2-4.6454090319E-7*T^3+8.9042735735E-10*T^4-9.0790692686E-13*T^5+3.8457331488E-16*T^6'; '413.15' '553.75' '0.00401235783-2.10746715E-5*T^1+3.85772275E-8*T^2-2.39730284E-11*T^3'});
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('fununit', 'Pa*s');
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('pieces', {'273.15' '553.75' '12010.1471-80.4072879*T^1+0.309866854*T^2-5.38186884E-4*T^3+3.62536437E-7*T^4'});
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('fununit', 'J/(kg*K)');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('smooth', 'contd1');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('pieces', {'273.15' '293.15' '0.000063092789034*T^3-0.060367639882855*T^2+18.9229382407066*T-950.704055329848'; '293.15' '373.15' '0.000010335053319*T^3-0.013395065634452*T^2+4.969288832655160*T+432.257114008512'});
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('fununit', 'kg/m^3');
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('pieces', {'273.15' '553.75' '-0.869083936+0.00894880345*T^1-1.58366345E-5*T^2+7.97543259E-9*T^3'});
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('fununit', 'W/(m*K)');
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('table', {'273' '1403';  ...
'278' '1427';  ...
'283' '1447';  ...
'293' '1481';  ...
'303' '1507';  ...
'313' '1526';  ...
'323' '1541';  ...
'333' '1552';  ...
'343' '1555';  ...
'353' '1555';  ...
'363' '1550';  ...
'373' '1543'});
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('interp', 'piecewisecubic');
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('fununit', 'm/s');
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('funcname', 'alpha_p');
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('expr', '-1/rho(T)*d(rho(T),T)');
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('args', {'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('fununit', '1/K');
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('plotargs', {'T' '273.15' '373.15'});
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('funcname', 'gamma_w');
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('expr', '1+(T/Cp(T))*(alpha_p(T)*cs(T))^2');
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('args', {'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('fununit', '1');
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('plotargs', {'T' '273.15' '373.15'});
model.component('comp1').material('mat1').propertyGroup('def').func('an3').set('funcname', 'muB');
model.component('comp1').material('mat1').propertyGroup('def').func('an3').set('expr', '2.79*eta(T)');
model.component('comp1').material('mat1').propertyGroup('def').func('an3').set('args', {'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('an3').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('def').func('an3').set('fununit', 'Pa*s');
model.component('comp1').material('mat1').propertyGroup('def').func('an3').set('plotargs', {'T' '273.15' '553.75'});
model.component('comp1').material('mat1').propertyGroup('def').set('thermalexpansioncoefficient', '');
model.component('comp1').material('mat1').propertyGroup('def').set('bulkviscosity', '');
model.component('comp1').material('mat1').propertyGroup('def').set('thermalexpansioncoefficient', {'alpha_p(T)' '0' '0' '0' 'alpha_p(T)' '0' '0' '0' 'alpha_p(T)'});
model.component('comp1').material('mat1').propertyGroup('def').set('bulkviscosity', 'muB(T)');
model.component('comp1').material('mat1').propertyGroup('def').descr('thermalexpansioncoefficient_symmetry', '');
model.component('comp1').material('mat1').propertyGroup('def').descr('bulkviscosity_symmetry', '');
model.component('comp1').material('mat1').propertyGroup('def').set('dynamicviscosity', 'eta(T)');
model.component('comp1').material('mat1').propertyGroup('def').descr('dynamicviscosity_symmetry', '');
model.component('comp1').material('mat1').propertyGroup('def').set('ratioofspecificheat', 'gamma_w(T)');
model.component('comp1').material('mat1').propertyGroup('def').descr('ratioofspecificheat_symmetry', '');
model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', {'5.5e-6[S/m]' '0' '0' '0' '5.5e-6[S/m]' '0' '0' '0' '5.5e-6[S/m]'});
model.component('comp1').material('mat1').propertyGroup('def').descr('electricconductivity_symmetry', '');
model.component('comp1').material('mat1').propertyGroup('def').set('heatcapacity', 'Cp(T)');
model.component('comp1').material('mat1').propertyGroup('def').descr('heatcapacity_symmetry', '');
model.component('comp1').material('mat1').propertyGroup('def').set('density', 'rho(T)');
model.component('comp1').material('mat1').propertyGroup('def').descr('density_symmetry', '');
model.component('comp1').material('mat1').propertyGroup('def').set('thermalconductivity', {'k(T)' '0' '0' '0' 'k(T)' '0' '0' '0' 'k(T)'});
model.component('comp1').material('mat1').propertyGroup('def').descr('thermalconductivity_symmetry', '');
model.component('comp1').material('mat1').propertyGroup('def').set('soundspeed', 'cs(T)');
model.component('comp1').material('mat1').propertyGroup('def').descr('soundspeed_symmetry', '');
model.component('comp1').material('mat1').propertyGroup('def').addInput('temperature');

model.component('comp1').common('rot1').set('rotationType', 'rotationalVelocity');
model.component('comp1').common('rot1').set('rotationalVelocityExpression', 'generalRevolutionsPerTime');
model.component('comp1').common('rot1').set('revolutionsPerTime', 1);

model.component('comp1').physics('spf').feature('vf1').set('F', {'0'; '0'; '-mat1.def.rho*9.81[m/s^2]'});

model.component('comp1').mesh('mesh1').feature('size').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('size').set('table', 'cfd');
model.component('comp1').mesh('mesh1').feature('size').set('hmax', '1[mm]');
model.component('comp1').mesh('mesh1').feature('size').set('hmin', '1[mm]');
model.component('comp1').mesh('mesh1').feature('size').set('hnarrow', 0.7);
model.component('comp1').mesh('mesh1').feature('size').set('hgrad', 1.15);
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hmax', '1[mm]');
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hmin', '1[mm]');
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').run;
model.component('comp1').mesh('mesh2').contribute('physics/spf', false);

model.study.create('std1');
model.study('std1').create('frrot', 'FrozenRotor');
model.study('std1').create('sfs', 'StationaryFreeSurface');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').create('su1', 'StoreSolution');
model.sol('sol1').create('st2', 'StudyStep');
model.sol('sol1').create('v2', 'Variables');
model.sol('sol1').create('s2', 'Stationary');
model.sol('sol1').feature('s1').create('se1', 'Segregated');
model.sol('sol1').feature('s1').create('i1', 'Iterative');
model.sol('sol1').feature('s1').create('i2', 'Iterative');
model.sol('sol1').feature('s1').create('d1', 'Direct');
model.sol('sol1').feature('s1').create('d2', 'Direct');
model.sol('sol1').feature('s1').feature('se1').create('ss1', 'SegregatedStep');
model.sol('sol1').feature('s1').feature('se1').create('ss2', 'SegregatedStep');
model.sol('sol1').feature('s1').feature('se1').create('ll1', 'LowerLimit');
model.sol('sol1').feature('s1').feature('se1').feature.remove('ssDef');
model.sol('sol1').feature('s1').feature('i1').create('mg1', 'Multigrid');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').create('sc1', 'SCGS');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').create('sc1', 'SCGS');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol1').feature('s1').feature('i2').create('mg1', 'Multigrid');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('pr').create('sl1', 'SORLine');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('po').create('sl1', 'SORLine');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol1').feature('s1').feature.remove('fcDef');
model.sol('sol1').feature('s2').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s2').create('i1', 'Iterative');
model.sol('sol1').feature('s2').create('d1', 'Direct');
model.sol('sol1').feature('s2').feature('i1').create('mg1', 'Multigrid');
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('pr').create('sc1', 'SCGS');
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('po').create('sc1', 'SCGS');
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol1').feature('s2').feature.remove('fcDef');

model.result.dataset.create('surf1', 'Surface');
model.result.dataset('surf1').selection.set([1 2 3 14 17]);
model.result.create('pg1', 'PlotGroup3D');
model.result.create('pg2', 'PlotGroup3D');
model.result.create('pg3', 'PlotGroup3D');
model.result.create('pg4', 'PlotGroup3D');
model.result.create('pg6', 'PlotGroup3D');
model.result.create('pg8', 'PlotGroup3D');
model.result('pg1').create('slc1', 'Slice');
model.result('pg1').create('arwv1', 'ArrowVolume');
model.result('pg2').create('surf1', 'Surface');
model.result('pg2').create('con1', 'Contour');
model.result('pg2').feature('surf1').set('expr', '1');
model.result('pg2').feature('con1').set('expr', 'p');
model.result('pg3').create('surf1', 'Surface');
model.result('pg3').feature('surf1').set('expr', 'spf.Delta_wPlus');
model.result('pg4').create('arwv1', 'ArrowVolume');
model.result('pg4').create('slc1', 'Slice');
model.result('pg6').create('slc1', 'Slice');
model.result('pg6').feature('slc1').set('expr', '((mat1.def.eta(273))^3/mat1.def.rho(273)^3/ep)^0.25');
model.result('pg8').create('slc1', 'Slice');
model.result('pg8').feature('slc1').set('expr', 'spf.sr*mat1.def.eta(273)');
model.result.export.create('data1', 'Data');
model.result.export.create('data2', 'Data');
model.result.export.create('data3', 'Data');

model.study('std1').feature('frrot').set('mesh', {'geom1' 'mesh1'});
model.study('std1').feature('frrot').set('ngen', 5);
model.study('std1').feature('sfs').set('mesh', {'geom1' 'mesh1'});

model.sol('sol1').attach('std1');
model.sol('sol1').feature('s1').feature('aDef').set('cachepattern', true);
model.sol('sol1').feature('s1').feature('se1').set('maxsegiter', 400);
model.sol('sol1').feature('s1').feature('se1').set('segstabacc', 'segcflcmp');
model.sol('sol1').feature('s1').feature('se1').set('subinitcfl', 3);
model.sol('sol1').feature('s1').feature('se1').feature('ss1').label('Velocity u, Pressure p');
model.sol('sol1').feature('s1').feature('se1').feature('ss1').set('segvar', {'comp1_u' 'comp1_p'});
model.sol('sol1').feature('s1').feature('se1').feature('ss1').set('linsolver', 'i1');
model.sol('sol1').feature('s1').feature('se1').feature('ss1').set('subdamp', 0.5);
model.sol('sol1').feature('s1').feature('se1').feature('ss2').label('Turbulence variables');
model.sol('sol1').feature('s1').feature('se1').feature('ss2').set('segvar', {'comp1_k' 'comp1_ep'});
model.sol('sol1').feature('s1').feature('se1').feature('ss2').set('linsolver', 'i2');
model.sol('sol1').feature('s1').feature('se1').feature('ss2').set('subdamp', 0.35);
model.sol('sol1').feature('s1').feature('se1').feature('ss2').set('subtermconst', 'itertol');
model.sol('sol1').feature('s1').feature('se1').feature('ss2').set('subiter', 3);
model.sol('sol1').feature('s1').feature('se1').feature('ss2').set('subntolfact', 1);
model.sol('sol1').feature('s1').feature('se1').feature('ll1').set('lowerlimit', 'comp1.k 0 comp1.ep 0 ');
model.sol('sol1').feature('s1').feature('i1').label('GMG, fluid flow variables (spf)');
model.sol('sol1').feature('s1').feature('i1').set('nlinnormuse', true);
model.sol('sol1').feature('s1').feature('i1').set('maxlinit', 200);
model.sol('sol1').feature('s1').feature('i1').set('rhob', 20);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sc1').set('linesweeptype', 'ssor');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sc1').set('iter', 0);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sc1').set('scgsvertexrelax', 0.7);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sc1').set('linesweeptype', 'ssor');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sc1').set('iter', 1);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sc1').set('scgsvertexrelax', 0.7);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').feature('d1').set('pivotperturb', 1.0E-13);
model.sol('sol1').feature('s1').feature('i2').label('AMG, turbulence variables (spf)');
model.sol('sol1').feature('s1').feature('i2').set('nlinnormuse', true);
model.sol('sol1').feature('s1').feature('i2').set('maxlinit', 200);
model.sol('sol1').feature('s1').feature('i2').set('rhob', 20);
model.sol('sol1').feature('s1').feature('i2').feature('mg1').set('prefun', 'saamg');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').set('maxcoarsedof', 50000);
model.sol('sol1').feature('s1').feature('i2').feature('mg1').set('saamgcompwise', true);
model.sol('sol1').feature('s1').feature('i2').feature('mg1').set('usesmooth', false);
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('pr').feature('sl1').set('linesweeptype', 'ssor');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('pr').feature('sl1').set('iter', 0);
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('pr').feature('sl1').set('linerelax', 0.7);
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('pr').feature('sl1').set('linemethod', 'uncoupled');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('pr').feature('sl1').set('relax', 0.5);
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('po').feature('sl1').set('linesweeptype', 'ssor');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('po').feature('sl1').set('iter', 1);
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('po').feature('sl1').set('linerelax', 0.7);
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('po').feature('sl1').set('linemethod', 'uncoupled');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('po').feature('sl1').set('relax', 0.5);
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('cs').feature('d1').set('pivotperturb', 1.0E-13);
model.sol('sol1').feature('s1').feature('d1').label('Direct, fluid flow variables (spf)');
model.sol('sol1').feature('s1').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('s1').feature('d1').set('pivotperturb', 1.0E-13);
model.sol('sol1').feature('s1').feature('d2').label('Direct, turbulence variables (spf)');
model.sol('sol1').feature('s1').feature('d2').set('linsolver', 'pardiso');
model.sol('sol1').feature('s1').feature('d2').set('pivotperturb', 1.0E-13);
model.sol('sol1').feature('st2').set('studystep', 'sfs');
model.sol('sol1').feature('v2').set('initmethod', 'sol');
model.sol('sol1').feature('v2').set('initsol', 'sol1');
model.sol('sol1').feature('v2').set('solnum', 'auto');
model.sol('sol1').feature('v2').set('notsolmethod', 'sol');
model.sol('sol1').feature('v2').set('notsol', 'sol1');
model.sol('sol1').feature('v2').set('notsolnum', 'auto');
model.sol('sol1').feature('s2').feature('aDef').set('cachepattern', true);
model.sol('sol1').feature('s2').feature('fc1').set('linsolver', 'i1');
model.sol('sol1').feature('s2').feature('fc1').set('dtech', 'const');
model.sol('sol1').feature('s2').feature('fc1').set('maxiter', 300);
model.sol('sol1').feature('s2').feature('fc1').set('damp', 0.8);
model.sol('sol1').feature('s2').feature('fc1').set('stabacc', 'cflcmp');
model.sol('sol1').feature('s2').feature('fc1').set('initcfl', 3);
model.sol('sol1').feature('s2').feature('i1').label('GMG, fluid flow variables (spf)');
model.sol('sol1').feature('s2').feature('i1').set('nlinnormuse', true);
model.sol('sol1').feature('s2').feature('i1').set('maxlinit', 200);
model.sol('sol1').feature('s2').feature('i1').set('rhob', 20);
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('pr').feature('sc1').set('linesweeptype', 'ssor');
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('pr').feature('sc1').set('iter', 0);
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('pr').feature('sc1').set('scgsvertexrelax', 0.7);
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('po').feature('sc1').set('linesweeptype', 'ssor');
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('po').feature('sc1').set('iter', 1);
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('po').feature('sc1').set('scgsvertexrelax', 0.7);
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('cs').feature('d1').set('pivotperturb', 1.0E-13);
model.sol('sol1').feature('s2').feature('d1').label('Direct, fluid flow variables (spf)');
model.sol('sol1').feature('s2').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('s2').feature('d1').set('pivotperturb', 1.0E-13);
model.sol('sol1').runAll;

model.result.dataset('surf1').label('Exterior Walls');
model.result('pg1').label('Velocity (spf)');
model.result('pg1').set('frametype', 'spatial');
model.result('pg1').set('showlegendsmaxmin', true);
model.result('pg1').set('showlegendsunit', true);
model.result('pg1').feature('slc1').label('Slice');
model.result('pg1').feature('slc1').set('quickplane', 'xy');
model.result('pg1').feature('slc1').set('smooth', 'internal');
model.result('pg1').feature('slc1').set('resolution', 'normal');
model.result('pg1').feature('arwv1').set('scale', 0.0725446997392873);
model.result('pg1').feature('arwv1').set('scaleactive', false);
model.result('pg2').label('Pressure (spf)');
model.result('pg2').set('data', 'surf1');
model.result('pg2').set('frametype', 'spatial');
model.result('pg2').feature('surf1').label('Surface');
model.result('pg2').feature('surf1').set('titletype', 'none');
model.result('pg2').feature('surf1').set('coloring', 'uniform');
model.result('pg2').feature('surf1').set('color', 'gray');
model.result('pg2').feature('surf1').set('smooth', 'internal');
model.result('pg2').feature('surf1').set('resolution', 'normal');
model.result('pg2').feature('con1').label('Pressure');
model.result('pg2').feature('con1').set('number', 40);
model.result('pg2').feature('con1').set('levelrounding', false);
model.result('pg2').feature('con1').set('smooth', 'internal');
model.result('pg2').feature('con1').set('resolution', 'normal');
model.result('pg3').label('Wall Resolution (spf)');
model.result('pg3').set('data', 'surf1');
model.result('pg3').set('frametype', 'spatial');
model.result('pg3').feature('surf1').label('Wall Resolution');
model.result('pg3').feature('surf1').set('smooth', 'internal');
model.result('pg3').feature('surf1').set('resolution', 'normal');
model.result('pg4').label('Velocity (spf) 1');
model.result('pg4').set('frametype', 'spatial');
model.result('pg4').set('showlegendsmaxmin', true);
model.result('pg4').set('showlegendsunit', true);
model.result('pg4').feature('arwv1').set('scale', 0.0725446997392873);
model.result('pg4').feature('arwv1').set('scaleactive', false);
model.result('pg4').feature('slc1').set('quickplane', 'xy');
model.result('pg4').feature('slc1').set('quickznumber', 1);
model.result('pg4').feature('slc1').set('resolution', 'normal');
model.result('pg6').label('Kolmogorov length');
model.result('pg6').set('frametype', 'spatial');
model.result('pg6').set('showlegendsmaxmin', true);
model.result('pg6').set('showlegendsunit', true);
model.result('pg6').feature('slc1').label('Slice');
model.result('pg6').feature('slc1').set('descractive', true);
model.result('pg6').feature('slc1').set('descr', 'Kolmogorov length');
model.result('pg6').feature('slc1').set('quickplane', 'zx');
model.result('pg6').feature('slc1').set('quickynumber', 1);
model.result('pg6').feature('slc1').set('smooth', 'internal');
model.result('pg6').feature('slc1').set('resolution', 'normal');
model.result('pg8').label('Shear stress');
model.result('pg8').set('frametype', 'spatial');
model.result('pg8').set('showlegendsmaxmin', true);
model.result('pg8').set('showlegendsunit', true);
model.result('pg8').feature('slc1').label('Slice');
model.result('pg8').feature('slc1').set('descractive', true);
model.result('pg8').feature('slc1').set('descr', 'Shear stress');
model.result('pg8').feature('slc1').set('quickplane', 'zx');
model.result('pg8').feature('slc1').set('quickynumber', 1);
model.result('pg8').feature('slc1').set('smooth', 'internal');
model.result('pg8').feature('slc1').set('resolution', 'normal');
model.result.export('data1').label('Data 1 -velocity');
model.result.export('data1').set('expr', {'u' 'v' 'w'});
model.result.export('data1').set('unit', {'m/s' 'm/s' 'm/s'});
model.result.export('data1').set('descr', {'Velocity field, x component' 'Velocity field, y component' 'Velocity field, z component'});
model.result.export('data1').set('filename', 'D:\Jiajia Chen\Cell-based Meat\CFD export for ABM\Tetrahedra mesh\Velocity.txt');
model.result.export('data2').label('Data 2 - Kolmogorov length');
model.result.export('data2').set('expr', {'((mat1.def.eta(273))^3/mat1.def.rho(273)^3/ep)^0.25'});
model.result.export('data2').set('unit', {'m'});
model.result.export('data2').set('descr', {''});
model.result.export('data2').set('filename', 'D:\Jiajia Chen\Cell-based Meat\CFD export for ABM\Tetrahedra mesh\Kolmogorov length - 60 rpm.txt');
model.result.export('data3').label('Data 2 - Shear stress');
model.result.export('data3').set('expr', {'spf.sr*mat1.def.eta(273)'});
model.result.export('data3').set('unit', {'Pa'});
model.result.export('data3').set('descr', {''});
model.result.export('data3').set('filename', 'D:\Jiajia Chen\Cell-based Meat\CFD export for ABM\Tetrahedra mesh\Shear Stress - 60 rpm.txt');

out = model;
