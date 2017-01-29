within ;
package EnergyStation
  package BoreholeField

    function erf "error function "
      input Real u;
      output Real y;
    external"C" y = erf(u) annotation (Include="#include <erf.c>",
                       IncludeDirectory="modelica://EnergyStation/C-Sources");
      annotation (Documentation(info="<html>
<p>
Error function that recieves the input value <code>u</code> ,
and that returns the value <code>y</code> by calling the external function.
</p>
</html>",     revisions="<html>
<ul><li>
July 3, 2014, by Liang Yin:<br/>
First implementation.
</li>
</ul>
</html>"));

    end erf;

    model Testerf
       Real x=time;
       Real y;
       Real iy;
       Real diy;
       Real Iis;
       Real Eis;
       Real sEis;
      // Real a[2,2]={{1,2},{3,4}};
      // Real b[1,2]={{5,6}};
      // Real c[:,2]=cat(1,a,b);

     //Real a[:,2];

      //    t[1]=Ei((1/(1e-6*(1e-3))),positions={{0,0},{5,0},{5,5},{0,5}});
      //    t[1]=Ei((1/(1e-6*(1))),positions={{0,0},{5,0},{5,5},{0,5}});
     // a:=cat(1,a,{{time,Eis}});
     //   paramter Real coeff[6]=BoreholeField.coeffEi();
      parameter Real coeff[:]=BoreholeField.coeffEi(r_b=0.12,positions={{0,0},{5,0},{5,5},{0,5}});
      parameter Real ek=BoreholeField.TestcoeffEiFunc(2e7, coeff);
      Real Eif;
    equation
      Eif = BoreholeField.TestcoeffEiFunc((time + 1e-3), coeff);
       y=erf(x);
       iy=ierf(x);
       der(diy)=y;
       Iis=I_ls(x);
       //   Eis=Ei((time+0.0001),positions={{0,0},{500,0}});
        Eis=Ei((1/(1e-6*(time+1e-3))),r_b=0.12,positions={{0,0},{5,0},{5,5},{0,5}});
       sEis=Buildings.Fluid.HeatExchangers.Boreholes.BaseClasses.powerSeries(1/(1e-6*(time+20000))*0.1^2,3);
    end Testerf;

    model TestEi "Compare the original model with my model of the g-function"
       Real x=time;
       parameter Real alpha=1e-6;

       Real mEis;
       Real sEis;

    //  parameter Real coeff[:]=BoreholeField.coeffEi(H=100,r_b=0.12,positions={{0,0},{5,0},{5,5},{0,5}});
    parameter Real coeff[:]=BoreholeField.coeffEi(H=1000,r_b=0.12,positions={{0,0}});
      //BoreholeField.coeffEiFunc cf(lookTable(table=coeff));

    //  parameter Real ek=BoreholeField.coeffEiFunc(2e7,coeff);
    parameter Real ekp=Modelica.Math.Vectors.interpolate(BoreholeField.Constants.timePoint,coeff,(1e10)*alpha/
          BoreholeField.Constants.coeffEiP);
    parameter Real ek=Modelica.Math.Vectors.interpolate(BoreholeField.Constants.timePoint,coeff,(1e9*10)*alpha/
          BoreholeField.Constants.coeffEiP);
    parameter Real ek2=Ei((1/(alpha*((1e9)*10))),H=1000,r_b=0.12,positions={{0,0}});
     parameter Real ek3=Ei((1/(alpha*((1e10)*1))),H=1000,r_b=0.12,positions={{0,0}});
       Real Eif;
    equation
       //cf.u=((time+1e-3)*4*alpha/BoreholeField.Constants.coeffEiP);
       //Eif = cf.y;
    //   Eif = Modelica.Math.Vectors.interpolate(BoreholeField.Constants.timePoint,coeff,((time+1e-3)*4*alpha/BoreholeField.Constants.coeffEiP));
       // Eif=cf(((time+1e-3)*4*alpha/BoreholeField.Constants.coeffEiP),   coeff);
      Eif = Modelica.Math.Vectors.interpolate(
            BoreholeField.Constants.timePoint,
            coeff,
            ((time + 1e4)*4*alpha/BoreholeField.Constants.coeffEiP));
     //   mEis=Ei((1/(4*alpha*(time+1e-3))),H=1000,r_b=0.12,positions={{0,0},{5,0},{5,5},{0,5}});
       mEis=Ei((1/(4*alpha*(time+1e4))),H=1000,r_b=0.12,positions={{0,0}});
       sEis=Buildings.Fluid.HeatExchangers.Boreholes.BaseClasses.powerSeries(1/(4*1e-6*(time+1e4))*0.12^2,2);
         annotation (Documentation(info="<html>
<p>
The long time error come from the model object difference that the original is infinite whereas my model is finite height
</p>
</html>",
    revisions="<html>
<ul><li>
July 3, 2014, by Liang Yin:<br/>
First implementation.
</li>
</ul>
</html>"));
    end TestEi;

    function ierf
      input Real u;
      output Real y;
    algorithm
      y:=u*erf(u)-(1-exp(-u^2))/sqrt(Modelica.Constants.pi);
      annotation (Documentation(info="<html>
<p>
Integral error function that recieves the input value <code>u</code> ,
and that returns the integral error function value from <code>0</code> to <code>u</code> .
</p>
</html>",     revisions="<html>
<ul><li>
July 3, 2014, by Liang Yin:<br/>
First implementation.
</li>
</ul>
</html>"));
    end ierf;

    function I_ls
      input Real u;
      output Real y;
    algorithm
      y:=4*ierf(u)-ierf(2*u);
        annotation (Documentation(info="<html>
<p>
The function that recieves the input value <code>u</code> ,
and that returns the value <code>y</code> .
</p>
<h4>References</h4>
<dl>
 <dt>Damien Picard, Lieve Helsen:</dt>
 <dd><b>Advanced Hybrid Model for Borefield Heat Exchanger Performance Evaluation, an Implementation in Modelica</b>.
 In Proceedings of the 10th International Modelica Conference, volume 1, pages 857-866, 2014.</dd>

</dl>
 
</html>",     revisions="<html>
<ul><li>
July 3, 2014, by Liang Yin:<br/>
First implementation.
</li>
</ul>
</html>"));
    end I_ls;

    function coeffEi
       input Real r_b;
       input Real positions[:,:];
       input Real H;
       output Real xt[N];
     //  import Modelica.Utilities.Streams.print;
       constant Integer N=BoreholeField.Constants.N;

       //Real yt[N];
       //Real matrixt[N,N];
       //Real tem[N];

    algorithm
      for i in 1:N loop
     //   yt[i]:=Ei((1/(1e-6*(base^(i)))),positions={{0,0},{5,0},{5,5},{0,5}});
        xt[i] := Ei(
              (1/(BoreholeField.Constants.coeffEiP*BoreholeField.Constants.timePoint[
            i])),
              H=H,
              r_b=r_b,
              positions=positions);
      //  print("yt "+String(i)+" "+String(yt[i]));
      //  xt[i,1]:= log(BoreholeField.Constants.timePoint[i]);
      //  xt[i,2]:=yt[i];
      end for;

    /*
  for i in 1:size(yt,1) loop
    for j in 1:size(yt,1) loop
      matrixt[i,j]:=(log(BoreholeField.Constants.timePoint[i]))^(j - 1);
  end for;
  end for;

//  xt:=Modelica.Math.Matrices.solve(matrixt,yt);

xt:=Modelica.Math.Matrices.leastSquares(matrixt,yt);
*/
    //tem:=matrixt*xt;
    /*
    for i in 1:size(yt,1) loop
    for j in 1:size(yt,1) loop
      tem[i]:=tem[i]+(log(timePoint[i]))^(j-1)*xt[j];
    end for;
    end for;
  for i in 1:size(yt,1) loop
    print("xt "+String(i)+" "+String(xt[i]));
    print("tem "+String(i)+" "+String(tem[i]));
    end for;
    */
    end coeffEi;

    function TestcoeffEiFunc
        input Real u;
        output Real y;
       // Modelica.Blocks.Tables.CombiTable1Ds lookTable(table=fill(0,0,2));

       Modelica.Math.Vectors.interpolate gg;
    algorithm
      /*
     y:=0;
    for j in 1:size(cE,1) loop
      y:=y+(log(u))^(j-1)*cE[j];
      end for;
      */
          lookTable.u:=log(u);
          y:=lookTable.y[1];
    end TestcoeffEiFunc;

    function Ei "Ei used to compute far-field temperature"
      input Real u "u=1/(4*alpha*t)";
      input Real r_b = 0.1 "The radium of the borehole";
      input Real H=100;
      input Real[:,2] positions "positions of the borehole";
      output Real W "The integal value";

    protected
      parameter Real tolerance= 1e-9;
    //  parameter Real H=150;
      parameter Real inf=1e7;

      Real s;

    algorithm
      s:=sqrt(u);

      W:= Modelica.Math.Nonlinear.quadratureLobatto(
          function integalFuncEi(H=H,r_b=r_b,positions=positions),
          s,
          inf,
          tolerance);

    end Ei;

    function integalFuncEi "The intergal term for Ei"
          extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
          input Real[:,2] positions;
          input Real r_b;
          input Real H;
          import Modelica.Utilities.Streams.print;
    protected
        Real sumRS;
        Real r_ij;
        Real multi;

    algorithm
        sumRS:=0;
        multi:=2;
        for i in 1:size(positions,1) loop
           for j in 1:i loop
              if i==j then
                r_ij:=r_b;
                multi:=1;
              else
                r_ij:=sqrt((positions[i,1]-positions[j,1])^2+(positions[i,2]-positions[j,2])^2);
                multi:=2;
              end if;

              sumRS:=sumRS+ multi*exp(-r_ij^2*u^2);
            end for;
        end for;

         //          print("sumRS1 = "+String(sumRS));
         //     print("sumRS = "+String(sumRS));
        y:=(sumRS/size(positions,1))*I_ls(H*u)/(H*u^2);
    end integalFuncEi;

    function integalFuncEiO "The intergal term for Ei"
          extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
          input Real H;
          input Real r_b;
          input Real[:,2] positions;
    protected
        Real sumRS;
        Real r_ij;
        Real multi;

    algorithm
        sumRS:=0;
        multi:=2;
        for i in 1:size(positions,1) loop
           for j in 1:size(positions,1) loop
      //     for j in 1:i loop
              if i==j then
                r_ij:=r_b;
                multi:=1;
              else
                r_ij:=sqrt((positions[i,1]-positions[j,1])^2+(positions[i,2]-positions[j,2])^2);
                multi:=1;
              end if;
              sumRS:=sumRS+ multi*exp(-r_ij^2*u^2);

            end for;
        end for;

        y:=(sumRS/size(positions,1))*I_ls(H*u)/(H*u^2);
    end integalFuncEiO;

    function temperatureDrop
      "Calculate the temperature drop of the soil at the external boundary of the cylinder"
    input Buildings.Fluid.HeatExchangers.Boreholes.BaseClasses.ExtendableArray table
        "External object that contains the history terms of the heat flux";
    input Integer iSam(min=1)
        "Counter for how many time the model was sampled. Defined as iSam=1 when called at t=0";

    input Modelica.SIunits.HeatFlowRate Q_flow
        "Heat flow rate to be stored in the external object";
    input Modelica.SIunits.Time samplePeriod "Period between two samples";
    input Modelica.SIunits.Radius rExt "External radius of the cylinder";
    input Modelica.SIunits.Height hSeg "Height of the cylinder";
    input Modelica.SIunits.ThermalConductivity k
        "Thermal conductivity of the soil";
    input Modelica.SIunits.Density d "Density of the soil";
    input Modelica.SIunits.SpecificHeatCapacity c
        "Specific heat capacity of the soil";
        input Real[:] coeff "Temperature value at discrete time";
    output Modelica.SIunits.TemperatureDifference dT
        "Temperature drop of the soil";
    protected
     Modelica.SIunits.Time minSamplePeriod= rExt^2/(4*(k/c/d)*3.8)
        "Minimal lenght of the sampling period";
     Modelica.SIunits.HeatFlowRate QL_flow
        "Intermediate variable for heat flow rate at the lower bound of the time interval";
     Modelica.SIunits.HeatFlowRate QU_flow
        "Intermediate variable for heat flow rate at the upper bound of the time interval";

    algorithm
      assert(rExt*rExt/(4*(k/c/d)*samplePeriod)<=3.8,
      "The samplePeriod has to be bigger than " + String(minSamplePeriod) + " for convergence purpose.
  samplePeriod = "    + String(samplePeriod));
      if iSam == 1 then
        // First call, at t=0
        dT := 0;
        QL_flow := Buildings.Fluid.HeatExchangers.Boreholes.BaseClasses.exchangeValues(
                           table=table, iX=iSam, x=Q_flow, iY=iSam);
      else
        dT := 0;
        // The first evaluation is at iSam=2, in which case we have one term of the sum,
        // and t=samplePeriod=(iSam-1)*samplePeriod
       for i in 1:(iSam-1) loop
          QL_flow := Buildings.Fluid.HeatExchangers.Boreholes.BaseClasses.exchangeValues(
            table=table, iX=iSam, x=Q_flow, iY=iSam+1-i);
          QU_flow := Buildings.Fluid.HeatExchangers.Boreholes.BaseClasses.exchangeValues(
            table=table, iX=iSam, x=Q_flow, iY=iSam-i);
          // The division by hSeg is because QU_flow and QL_flow are in [W], but the equation
          // requires [W/m], i.e., heat flow rate per unit length of the line source.
          dT := dT + 1/(4*Modelica.Constants.pi*k)*
            Modelica.Math.Vectors.interpolate(
                BoreholeField.Constants.timePoint,
                coeff,
                ((i*samplePeriod)*4*k/(c*d)/BoreholeField.Constants.coeffEiP))*
            (QL_flow - QU_flow)/hSeg;
        //           BoreholeField.Ei(u=c*d/(4*k*i*samplePeriod),r_b=rExt,positions={{0,0},{5,0},{0,5},{5,5}})*

       //     Buildings.Fluid.HeatExchangers.Boreholes.BaseClasses.powerSeries(
       //        u=c*d/(4*k*i*samplePeriod)*rExt^2, N=10)*
     end for;
      end if;
    annotation(Documentation(info="<html>
<p>
This function calculates the temperature drop of the soil at the outer boundary of the cylinder.
The analytical formula of Hart and Couvillion (1986) for constant 
heat extraction is adapted to a non-constant heat flux. 
To adapt the formula for a variable rate of heat extraction,
different constant heat extraction rates, starting at different time instances, 
are super-imposed.
To obtain the temperature drop at the time <i>t=n*&Delta;t</i>, the effects of 
constant rate of heat extractions are super-imposed as
</p>
<p align=\"center\" style=\"font-style:italic;\">
&Delta;T ( r , t=n &Delta;t )= 1 &frasl; ( 4 &pi; k ) 
  &sum;  W(u(r, t= i &Delta;t)) (q<sub>n-i+1</sub>-q<sub>n-i</sub>),
</p>
<p>
where <i>r</i> is the radius for which the temperature is computed,
<i>k</i> is the thermal conductivity of the material,
<i>W</i> is a solution of the heat conduction in polar coordinates and
<i>q<sub>i</sub>=Q<sub>i</sub>/h</i> is
the specific rate of heat extraction per unit lenght at time 
<i>t=i &Delta;t</i>.
The value of 
<i>W</i> is obtained using
</p>
<p align=\"center\" style=\"font-style:italic;\">
W(u)=[-0.5772 - ln(u) + u - u<sup>2</sup>/(2 &nbsp; 2!) +u<sup>3</sup>/(3 &nbsp; 3!) - u<sup>4</sup>/(4 &nbsp; 4!) + ....].
</p>
<p>
where
<i>u(r,t)= c &rho; r<sup>2</sup> &frasl; (4 t k) </i>,
<i>&rho;</i> is the mass density and
<i>c</i> is the specific heat capacity per unit mass.
</p>
<h4>Implementation</h4>
<p>
The rate of heat flow <i>Q<sub>i</sub></i> is obtained from the function
<a href=\"modelica://Buildings.Fluid.HeatExchangers.Boreholes.BaseClasses.exchangeValues\">
Buildings.Fluid.HeatExchangers.Boreholes.BaseClasses.exchangeValues</a>.
</p>
<h4>References</h4>
<p>
Hart and Couvillion, (1986). <i>Earth Coupled Heat Transfer.</i>
Publication of the National Water Well Association. 
</p>
</html>",
    revisions="<html>
<ul>
<li>
July 27, 2011, by Pierre Vigouroux:<br/>
First implementation.
</li>
</ul>
</html>"));
    end temperatureDrop;

    model SingleUTubeBoundaryCondition
      "Prescribed temperature at the outer boundary of a single U tube borehole"
      replaceable parameter Buildings.HeatTransfer.Data.Soil.Generic matSoi
        "Thermal properties of the soil"
         annotation (choicesAllMatching=true);
        parameter Real[:,2] positionsBor "Positions of the boreholes";
       //     parameter Real coeff[:]=BoreholeField.coeffEi(r_b=rExt,H=hBor,positions=positionsBor);
       parameter Real coeff[:]=(DataFiles.readCSVmatrix("gFunction.csv"))[:,1];
        parameter Modelica.SIunits.Height hBor "Total height of the borehole";

      parameter Modelica.SIunits.Radius rExt=3
        "Distance from the brine where the calculation is performed";
      parameter Modelica.SIunits.Height hSeg=10 "Height of the segment";
      parameter Modelica.SIunits.Temperature TExt_start=283.15
        "Initial external temperature";
      parameter Modelica.SIunits.Time samplePeriod=604800
        "Period between two samples";
      Buildings.Fluid.HeatExchangers.Boreholes.BaseClasses.ExtendableArray table=
          Buildings.Fluid.HeatExchangers.Boreholes.BaseClasses.ExtendableArray()
        "Extentable array, used to store history of rate of heat flows";
      Modelica.SIunits.HeatFlowRate QAve_flow
        "Average heat flux over a time period";

      Modelica.Blocks.Interfaces.RealInput Q_flow(unit="W")
        "Heat flow rate at the center of the borehole, positive if heat is added to soil"
         annotation (Placement(transformation(extent={{-120,-100},{-80,-60}})));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port "Heat port"
        annotation (Placement(transformation(extent={{86,-10},
                {106,10}},      rotation=0), iconTransformation(extent={{86,-10},
                {106,10}})));
    protected
      final parameter Modelica.SIunits.SpecificHeatCapacity c= matSoi.c
        "Specific heat capacity of the soil";
      final parameter Modelica.SIunits.ThermalConductivity k= matSoi.k
        "Thermal conductivity of the soil";
      final parameter Modelica.SIunits.Density d = matSoi.d
        "Density of the soil";
      Modelica.SIunits.Energy UOld "Internal energy at the previous period";
      Modelica.SIunits.Energy U
        "Current internal energy, defined as U=0 for t=tStart";
      final parameter Modelica.SIunits.Time startTime(fixed=false)
        "Start time of the simulation";
      Integer iSam(min=1)
        "Counter for how many time the model was sampled. Defined as iSam=1 when called at t=0";
    initial algorithm
      U         := 0;
      UOld      := 0;
      startTime := time;
      iSam      := 1;
    equation
      der(U) = Q_flow;
    algorithm
      when initial() or sample(startTime,samplePeriod) then
        QAve_flow := (U-UOld)/samplePeriod;
        UOld      := U;
        port.T := TExt_start + BoreholeField.temperatureDrop(
              table=table,
              iSam=iSam,
              Q_flow=QAve_flow,
              samplePeriod=samplePeriod,
              rExt=rExt,
              hSeg=hSeg,
              k=k,
              d=d,
              c=c,
              coeff=coeff);
        iSam := iSam+1;
      end when;

    annotation (
        Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                100,100}}), graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              pattern=LinePattern.None,
              fillColor={159,159,223},
              fillPattern=FillPattern.Backward),
            Line(
              points={{-102,0},{64,0}},
              color={191,0,0},
              thickness=0.5),
            Text(
              extent={{0,0},{-100,-100}},
              lineColor={0,0,0},
              textString="K"),
            Text(
              extent={{-150,150},{150,110}},
              textString="%name",
              lineColor={0,0,255}),
            Polygon(
              points={{40,-18},{40,22},{80,2},{40,-18}},
              lineColor={191,0,0},
              fillColor={191,0,0},
              fillPattern=FillPattern.Solid)}),
              Documentation(info="<html>
<p>
This model computes the temperature boundary condition at the outer boundary of the borehole.
It takes as an input the heat flow rate at the center of the borehole.
This heat flow rate is averaged over the sample period.
At each sampling interval, typically every one week, a new temperature boundary conditions is computed using
the analytical solution to a line source heat transfer problem.
</p>
<h4>Implementation</h4>
<p>
The computation of the temperature change of the boundary is computed using the function
<a href=\"modelica://Buildings.Fluid.HeatExchangers.Boreholes.BaseClasses.temperatureDrop\">
Buildings.Fluid.HeatExchangers.Boreholes.BaseClasses.temperatureDrop</a>.
</p>
</html>",
    revisions="<html>
<ul>
<li>
September 27, 2013, by Michael Wetter:<br/>
Moved assignment of <code>startTime</code> to <code>initial algorithm</code> section
to avoid an error in OpenModelica.
</li>
<li>
November 3 2011, by Michael Wetter:<br/>
Moved <code>der(U) := Q_flow;</code> from the algorithm section to the equation section
as this assignment does not conform to the Modelica specification.
</li>
<li>
September 9 2011, by Michael Wetter:<br/>
Moved <code>equation</code> section into <code>algorithm</code> section to make sure that the equations
in the <code>when</code> block are ordered correctly.
</li>
<li>
July 28 2011, by Pierre Vigouroux:<br/>
First implementation.
</li>
</ul>
</html>"));
    end SingleUTubeBoundaryCondition;

    model BoreholeSegment "Vertical segment of a borehole"
      extends Buildings.Fluid.Interfaces.PartialFourPortInterface(
         redeclare final package Medium1 = Medium,
         redeclare final package Medium2 = Medium,
         final m1_flow_nominal = m_flow_nominal,
         final m2_flow_nominal = m_flow_nominal,
         final m1_flow_small =   m_flow_small,
         final m2_flow_small =   m_flow_small,
         final allowFlowReversal1=allowFlowReversal,
         final allowFlowReversal2=allowFlowReversal);
      extends Buildings.Fluid.Interfaces.TwoPortFlowResistanceParameters;
      extends Buildings.Fluid.Interfaces.LumpedVolumeDeclarations(T_start=TFil_start);
      replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
        "Medium in the component" annotation (choicesAllMatching=true);
      replaceable parameter Buildings.HeatTransfer.Data.Soil.Generic matSoi
        "Thermal properties of soil"
        annotation (choicesAllMatching=true, Dialog(group="Soil"),
        Placement(transformation(extent={{2,70},{22,90}})));
      replaceable parameter
        Buildings.HeatTransfer.Data.BoreholeFillings.Generic                     matFil
        "Thermal properties of the filling material"
        annotation (choicesAllMatching=true, Dialog(group="Filling material"),
        Placement(transformation(extent={{-68,70},{-48,90}})));

      parameter Modelica.SIunits.MassFlowRate m_flow_nominal
        "Nominal mass flow rate"
        annotation(Dialog(group = "Nominal condition"));
      parameter Modelica.SIunits.MassFlowRate m_flow_small(min=0) = 1E-4*abs(m_flow_nominal)
        "Small mass flow rate for regularization of zero flow"
        annotation(Dialog(tab = "Advanced"));
      parameter Boolean homotopyInitialization = true
        "= true, use homotopy method"
        annotation(Evaluate=true, Dialog(tab="Advanced"));

      parameter Modelica.SIunits.Radius rTub=0.02 "Radius of the tubes"
        annotation (Dialog(group="Tubes"));
      parameter Modelica.SIunits.ThermalConductivity kTub=0.5
        "Thermal conductivity of the tubes" annotation (Dialog(group="Tubes"));
      parameter Modelica.SIunits.Length eTub=0.002 "Thickness of the tubes"
        annotation (Dialog(group="Tubes"));
      parameter Modelica.SIunits.Temperature TFil_start=283.15
        "Initial temperature of the filling material"
        annotation (Dialog(group="Filling material"));
        parameter Modelica.SIunits.Height hBor "Total height of the borehole";

      parameter Modelica.SIunits.Radius rExt=3
        "Radius of the soil used for the external boundary condition"
        annotation (Dialog(group="Soil"));
      parameter Modelica.SIunits.Temperature TExt_start=283.15
        "Initial far field temperature" annotation (Dialog(group="Soil"));
      parameter Integer nSta(min=1) = 10
        "Number of state variables in the soil"
        annotation (Dialog(group="Soil"));
      parameter Modelica.SIunits.Time samplePeriod=604800
        "Sample period for the external boundary condition"
        annotation (Dialog(group="Soil"));
        parameter Real[:,2] positionsBor "Positions of the boreholes";
      parameter Modelica.SIunits.Radius rBor=0.1 "Radius of the borehole";
      parameter Modelica.SIunits.Height hSeg "Height of the element";
      parameter Modelica.SIunits.Length xC=0.05
        "Shank spacing, defined as the distance between the center of a pipe and the center of the borehole";
    //  parameter Real B0=17.44 "Shape coefficient for grout resistance";
    //  parameter Real B1=-0.6052 "Shape coefficient for grout resistance";

     parameter Boolean allowFlowReversal = true
        "= true to allow flow reversal, false restricts to design direction (port_a -> port_b)"
        annotation(Dialog(tab="Assumptions"), Evaluate=true);

     Buildings.Fluid.HeatExchangers.Boreholes.BaseClasses.HexInternalElement pipFil(
        redeclare final package Medium = Medium,
        final matFil=matFil,
        final matSoi=matSoi,
        final hSeg=hSeg,
        final rTub=rTub,
        final eTub=eTub,
        final kTub=kTub,
        final kSoi=matSoi.k,
        final xC=xC,
        final rBor=rBor,
        final TFil_start=TFil_start,
        final m1_flow_nominal=m_flow_nominal,
        final m2_flow_nominal=m_flow_nominal,
        final dp1_nominal=dp_nominal,
        final dp2_nominal=0,
        final from_dp1=from_dp,
        final from_dp2=from_dp,
        final linearizeFlowResistance1=linearizeFlowResistance,
        final linearizeFlowResistance2=linearizeFlowResistance,
        final deltaM1=deltaM,
        final deltaM2=deltaM,
        final m1_flow_small=m_flow_small,
        final m2_flow_small=m_flow_small,
        final allowFlowReversal1=allowFlowReversal,
        final allowFlowReversal2=allowFlowReversal,
        final homotopyInitialization=homotopyInitialization,
        final energyDynamics=energyDynamics,
        final massDynamics=massDynamics,
        final p1_start=p_start,
        T1_start=T_start,
        X1_start=X_start,
        C1_start=C_start,
        C1_nominal=C_nominal,
        final p2_start=p_start,
        T2_start=T_start,
        X2_start=X_start,
        C2_start=C_start,
        C2_nominal=C_nominal)
        "Internal part of the borehole including the pipes and the filling material"
        annotation (Placement(transformation(extent={{-70,-10},{-50,10}})));
    //    final B0=B0,
    //    final B1=B1,
      Buildings.HeatTransfer.Conduction.SingleLayerCylinder soi(
        final material=matSoi,
        final h=hSeg,
        final nSta=nSta,
        final r_a=rBor,
        final r_b=rExt,
        final steadyStateInitial=false,
        final TInt_start=TFil_start,
        final TExt_start=TExt_start) "Heat conduction in the soil"
        annotation (Placement(transformation(extent={{0,-10},{20,10}})));
         BoreholeField.SingleUTubeBoundaryCondition TBouCon(
        final matSoi=matSoi,
        final positionsBor=positionsBor,
        final hBor=hBor,
        final rExt=rExt,
        final hSeg=hSeg,
        final TExt_start=TExt_start,
        final samplePeriod=samplePeriod)
        "Thermal boundary condition for the far-field"
        annotation (Placement(transformation(extent={{48,-10},{68,10}})));
    protected
      Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heaFlo
        annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));
    equation
      connect(pipFil.port_b1, port_b1)
                                    annotation (Line(
          points={{-50,6},{-40,6},{-40,60},{100,60}},
          color={0,127,255},
          smooth=Smooth.None));
      connect(pipFil.port_a2, port_a2)
                                    annotation (Line(
          points={{-50,-6},{-40,-6},{-40,-60},{100,-60}},
          color={0,127,255},
          smooth=Smooth.None));
      connect(pipFil.port_b2, port_b2)
                                    annotation (Line(
          points={{-70,-6},{-80,-6},{-80,-60},{-100,-60}},
          color={0,127,255},
          smooth=Smooth.None));
      connect(pipFil.port, heaFlo.port_a)
                                         annotation (Line(
          points={{-50,6.10623e-16},{-45,6.10623e-16},{-45,1.22125e-15},{-40,
              1.22125e-15},{-40,6.10623e-16},{-30,6.10623e-16}},
          color={191,0,0},
          smooth=Smooth.None));
      connect(heaFlo.port_b, soi.port_a) annotation (Line(
          points={{-10,6.10623e-16},{-7.5,6.10623e-16},{-7.5,1.22125e-15},{-5,
              1.22125e-15},{-5,6.10623e-16},{-5.55112e-16,6.10623e-16}},
          color={191,0,0},
          smooth=Smooth.None));
      connect(soi.port_b, TBouCon.port) annotation (Line(
          points={{20,6.10623e-16},{30,6.10623e-16},{30,20},{80,20},{80,6.10623e-16},
              {67.6,6.10623e-16}},
          color={191,0,0},
          smooth=Smooth.None));
      connect(port_a1, pipFil.port_a1) annotation (Line(
          points={{-100,60},{-80,60},{-80,6},{-70,6}},
          color={0,127,255},
          smooth=Smooth.None));
      connect(heaFlo.Q_flow, TBouCon.Q_flow) annotation (Line(
          points={{-20,-10},{-20,-20},{40,-20},{40,-8},{48,-8}},
          color={0,0,127},
          smooth=Smooth.None));
      annotation (
        Icon(graphics={
            Rectangle(
              extent={{-72,80},{68,-80}},
              lineColor={0,0,255},
              pattern=LinePattern.None,
              fillColor={95,95,95},
              fillPattern=FillPattern.Solid),
            Rectangle(
              extent={{88,54},{-88,64}},
              lineColor={0,0,255},
              pattern=LinePattern.None,
              fillColor={0,0,255},
              fillPattern=FillPattern.Solid),
            Rectangle(
              extent={{88,-64},{-88,-54}},
              lineColor={0,0,255},
              pattern=LinePattern.None,
              fillColor={0,0,255},
              fillPattern=FillPattern.Solid),
            Rectangle(
              extent={{-72,80},{68,68}},
              lineColor={0,0,0},
              fillColor={192,192,192},
              fillPattern=FillPattern.Backward),
            Rectangle(
              extent={{-72,-68},{68,-80}},
              lineColor={0,0,0},
              fillColor={192,192,192},
              fillPattern=FillPattern.Backward)}),
        Documentation(info="<html>
<p>
Horizontal layer that is used to model a U-tube borehole heat exchanger. 
This model combines three models, each simulating a different aspect 
of a borehole heat exchanger. 
</p>
<p>
The instance <code>pipFil</code> computes the heat transfer in the pipes and the filling material. 
This computation is done using the model
<a href=\"modelica://Buildings.Fluid.Boreholes.BaseClasses.HexInternalElement\">
Buildings.Fluid.Boreholes.BaseClasses.HexInternalElement</a>.
</p>
<p>
The instance <code>soi</code> computes transient and steady state heat transfer in the soil using a vertical cylinder.
The computation is done using the model <a href=\"modelica://Buildings.HeatTransfer.Conduction.SingleLayerCylinder\">
Buildings.HeatTransfer.Conduction.SingleLayerCylinder</a>.
</p>
<p>
The model <code>TBouCon</code> computes the far-field temperature boundary condition, i.e., the temperature at the outer
surface of the above cylindrical heat transfer computation.
The computation is done using the model
<a href=\"modelica://Buildings.Fluid.Boreholes.BaseClasses.TemperatureBoundaryCondition\">
Buildings.Fluid.Boreholes.BaseClasses.TemperatureBoundaryCondition</a>.
</p>
</html>",     revisions="<html>
<ul>
<li>
October 8, 2013, by Michael Wetter:<br/>
Removed parameter <code>show_V_flow</code>.
</li>
<li>
July 28 2011, by Pierre Vigouroux:<br/>
First implementation.
</li>
</ul>
</html>"));
    end BoreholeSegment;

    package Constants
        constant Integer N=55;
        constant Real coeffEiP= 1.0e-6;
        constant Real timePoint[:]={ 0.36,                3.6,               36,              360,              3*360,
                                     4*360,               6*360,             8*360,           9*360,            3600,
                                     3*3600,              6*3600,             8*3600,          12*3600,          16*3600,
                                     18*3600,             22*3600,           24*3600,         2*24*3600,        3*24*3600,
                                     5*24*3600,           6*24*3600,         18*24*3600,      60*24*3600,       100*24*3600,
                                     200*24*3600,         365*24*3600,       2*365*24*3600,   4*365*24*3600,    6*365*24*3600,
                                     10*365*24*3600,      14*365*24*3600,    20*365*24*3600,  40*365*24*3600,   60*365*24*3600,
                                     80*365*24*3600,      100*365*24*3600,   4.e9,            6.e9,             7.e9,
                                     8.e9,                1.e10,             2e10,            8.e10,           1.e11,
                                     2.e11,               3e11,              4e11,            5.8e11,           7e11,
                                     8.e11,               1.e12,             1.e13,           1.e14,            1.e15};

    end Constants;

    model twoPortWithParallel
      extends Modelica.Fluid.Interfaces.PartialTwoPort(
        port_a(p(start=Medium.p_default,
                 nominal=Medium.p_default)),
        port_b(p(start=Medium.p_default,
               nominal=Medium.p_default)));

       parameter Integer nParallel=1;

    equation
      port_a.p=port_b.p;

      port_a.m_flow+nParallel*port_b.m_flow=0;

      port_a.h_outflow=inStream(port_b.h_outflow);
      port_a.Xi_outflow=inStream(port_b.Xi_outflow);
      port_a.C_outflow=inStream(port_b.C_outflow);

      port_b.h_outflow=inStream(port_a.h_outflow);
      port_b.Xi_outflow=inStream(port_a.Xi_outflow);
      port_b.C_outflow=inStream(port_a.C_outflow);

    end twoPortWithParallel;

    model UTube "Single U-tube borehole heat exchanger"
      extends Buildings.Fluid.Interfaces.PartialTwoPortInterface(
        show_T=true);
      extends Buildings.Fluid.Interfaces.TwoPortFlowResistanceParameters(final
          computeFlowResistance=false, final linearizeFlowResistance=false);
      extends Buildings.Fluid.Interfaces.LumpedVolumeDeclarations;
      import Modelica.Constants;
         /*
    parameter Boolean positionsOnFile=false 
    "= true, if positions is defined on file"
        annotation (Dialog(group="BoreholeGroup"));
     
    parameter String fileName="NoName" "File where postions are stored"
    annotation (Dialog(
      group="BoreholeGroup",
      enable=positionsOnFile,
      loadSelector(filter="Text files (*.txt);;EXCEL CSV-files (*.csv)",
      caption="Open file in which positions is present")));
      */
     //   parameter Real[:,2] positionsBor= if positionsOnFile then DataFiles.readCSVmatrix("") else {{0,0}}
     parameter Real[:,2] positionsBor = DataFiles.readCSVmatrix("")
        "Positions of the boreholes"
        annotation (Dialog(group="BoreholeGroup"));
    //     parameter Real[:,2] positionsBor = DataFiles.readCSVmatrix("ma.csv")
        parameter Integer nBor = size(positionsBor,1) "Number of the boreholes"
         annotation (Dialog(group="BoreholeGroup"));
      replaceable parameter Buildings.HeatTransfer.Data.Soil.Generic matSoi
        "Thermal properties of soil"
        annotation (choicesAllMatching=true, Dialog(group="Soil"),
        Placement(transformation(extent={{2,70},{22,90}})));
      replaceable parameter
        Buildings.HeatTransfer.Data.BoreholeFillings.Generic                     matFil
        "Thermal properties of the filling material"
        annotation (choicesAllMatching=true, Dialog(group="Borehole"),
        Placement(transformation(extent={{-70,70},{-50,90}})));

      parameter Modelica.SIunits.Radius rTub=0.02 "Radius of the tubes"
        annotation(Dialog(group="Tubes"));
      parameter Modelica.SIunits.ThermalConductivity kTub=0.5
        "Thermal conductivity of the tube" annotation (Dialog(group="Tubes"));
      parameter Modelica.SIunits.Length eTub=0.002 "Thickness of a tube"
        annotation (Dialog(group="Tubes"));

      parameter Modelica.SIunits.Height hBor "Total height of the borehole"
        annotation(Dialog(group="Borehole"));
      parameter Integer nVer=10
        "Number of segments used for discretization in the vertical direction"
          annotation(Dialog(group="Borehole"));
      parameter Modelica.SIunits.Radius rBor=0.1 "Radius of the borehole";

      parameter Modelica.SIunits.Radius rExt=rBor+0.04
        "Radius of the soil used for the external boundary condition"
        annotation (Dialog(group="Soil"));
      parameter Integer nHor(min=1) = 10
        "Number of state variables in each horizontal layer of the soil"
        annotation (Dialog(group="Soil"));

      parameter Modelica.SIunits.Temperature TExt0_start=283.15
        "Initial far field temperature"
        annotation (Dialog(tab="Initial temperature", group="Soil"));
      parameter Modelica.SIunits.Temperature TExt_start[nVer]={if z[i] >= z0 then
          TExt0_start + (z[i] - z0)*dT_dz else TExt0_start for i in 1:nVer}
        "Temperature of the undisturbed ground"
        annotation (Dialog(tab="Initial temperature", group="Soil"));

      parameter Modelica.SIunits.Temperature TFil0_start=TExt0_start
        "Initial temperature of the filling material for h = 0...z0"
        annotation (Dialog(tab="Initial temperature", group="Filling material"));
      parameter Modelica.SIunits.Temperature TFil_start[nVer]=TExt_start
        "Temperature of the undisturbed ground"
        annotation (Dialog(tab="Initial temperature", group="Filling material"));

      parameter Modelica.SIunits.Height z0=10
        "Depth below which the temperature gradient starts"
        annotation (Dialog(tab="Initial temperature", group="Temperature profile"));
      parameter Real dT_dz(unit="K/m") = 0.01
        "Vertical temperature gradient of the undisturbed soil for h below z0"
        annotation (Dialog(tab="Initial temperature", group="Temperature profile"));

      parameter Modelica.SIunits.Time samplePeriod = 1e5
        "Sample period for the external boundary condition"
        annotation (Dialog(group="Soil"));
      parameter Modelica.SIunits.Length xC=0.05
        "Shank spacing, defined as the distance between the center of a pipe and the center of the borehole"
        annotation(Dialog(group="Borehole"));
    //  parameter Real B0=17.44 "Shape coefficient for grout resistance"
    //    annotation(Dialog(group="Borehole"));
    //  parameter Real B1=-0.605 "Shape coefficient for grout resistance"
    //    annotation(Dialog(group="Borehole"));
      parameter Boolean homotopyInitialization = true
        "= true, use homotopy method"
        annotation(Evaluate=true, Dialog(tab="Advanced"));
      BoreholeField.BoreholeSegment borHol[nVer](
        redeclare each final package Medium = Medium,
        each final matSoi=matSoi,
        each final matFil=matFil,
        each final hSeg=hBor/nVer,
        each final samplePeriod=samplePeriod,
        each final rTub=rTub,
        each final positionsBor=positionsBor,
        each final rBor=rBor,
        each final hBor=hBor,
        each final rExt=rExt,
        each final xC=xC,
        each final eTub=eTub,
        each final kTub=kTub,
        each final nSta=nHor,
        each final m_flow_nominal=m_flow_nominal,
        each final m_flow_small=m_flow_small,
        final dp_nominal={if i == 1 then dp_nominal else 0 for i in 1:nVer},
        TExt_start=TExt_start,
        TFil_start=TExt_start,
        each final homotopyInitialization=homotopyInitialization,
        each final show_T=show_T,
        each final computeFlowResistance=computeFlowResistance,
        each final from_dp=from_dp,
        each final linearizeFlowResistance=linearizeFlowResistance,
        each final deltaM=deltaM,
        each final energyDynamics=energyDynamics,
        each final massDynamics=massDynamics,
        each final p_start=p_start,
        each T_start=T_start,
        each X_start=X_start,
        each C_start=C_start,
        each C_nominal=C_nominal,
        each allowFlowReversal=allowFlowReversal)
        "Discretized borehole segments"
        annotation (Placement(transformation(extent={{-24,-10},{-4,10}})));
     //   each final B0=B0,
     //   each final B1=B1,

      Modelica.SIunits.Temperature Tdown[nVer] "Medium temperature in pipe 1";
      Modelica.SIunits.Temperature Tup[nVer] "Medium temperature in pipe 2";
    protected
      parameter Modelica.SIunits.Height z[nVer]={hBor/nVer*(i - 0.5) for i in 1:
          nVer} "Distance from the surface to the considered segment";
    public
      twoPortWithParallel twoPortWithParallel1(redeclare package Medium =  Medium,nParallel=nBor)
        annotation (Placement(transformation(extent={{-76,26},{-56,46}})));
      twoPortWithParallel twoPortWithParallel2(redeclare package Medium =  Medium,nParallel=nBor)
        annotation (Placement(transformation(extent={{58,28},{36,50}})));
    equation
      Tdown[:] = borHol[:].pipFil.vol1.heatPort.T;
      Tup[:] = borHol[:].pipFil.vol2.heatPort.T;
      connect(borHol[nVer].port_b1, borHol[nVer].port_a2) annotation (Line(
          points={{-4,6},{12,6},{12,-6},{-4,-6}},
          color={0,127,255},
          smooth=Smooth.None));
      for i in 1:nVer - 1 loop
        connect(borHol[i].port_b1, borHol[i + 1].port_a1) annotation (Line(
            points={{-4,6},{-4,20},{-24,20},{-24,6}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(borHol[i].port_a2, borHol[i + 1].port_b2) annotation (Line(
            points={{-4,-6},{-4,-20},{-24,-20},{-24,-6}},
            color={0,127,255},
            smooth=Smooth.None));
      end for;
      connect(port_b, port_b) annotation (Line(
          points={{100,0},{104,0},{104,0},{100,0}},
          color={0,127,255},
          smooth=Smooth.None));
      connect(port_b, twoPortWithParallel2.port_a) annotation (Line(
          points={{100,0},{90,0},{90,38},{58,38},{58,39}},
          color={0,127,255},
          smooth=Smooth.None));
      connect(twoPortWithParallel2.port_b, borHol[1].port_b2) annotation (Line(
          points={{36,39},{28,39},{28,-68},{-34,-68},{-34,-6},{-24,-6}},
          color={0,127,255},
          smooth=Smooth.None));
      connect(port_a, twoPortWithParallel1.port_a) annotation (Line(
          points={{-100,0},{-102,0},{-102,36},{-76,36}},
          color={0,127,255},
          smooth=Smooth.None));
      connect(twoPortWithParallel1.port_b, borHol[1].port_a1) annotation (Line(
          points={{-56,36},{-32,36},{-32,6},{-24,6}},
          color={0,127,255},
          smooth=Smooth.None));
        annotation (
        defaultComponentName="borehole",
        Icon(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-160,-100},{160,100}},
            grid={2,2},
            initialScale=0.5), graphics={
            Rectangle(
              extent={{-68,96},{72,-82}},
              lineColor={0,0,255},
              pattern=LinePattern.None,
              fillColor={95,95,95},
              fillPattern=FillPattern.Solid),
            Ellipse(
              extent={{-54,56},{-38,40}},
              lineColor={0,0,255},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid),
            Ellipse(
              extent={{-22,56},{-6,40}},
              lineColor={0,0,255},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid),
            Ellipse(
              extent={{12,56},{28,40}},
              lineColor={0,0,255},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid),
            Ellipse(
              extent={{42,56},{58,40}},
              lineColor={0,0,255},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid),
            Ellipse(
              extent={{-54,24},{-38,8}},
              lineColor={0,0,255},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid),
            Ellipse(
              extent={{-22,24},{-6,8}},
              lineColor={0,0,255},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid),
            Ellipse(
              extent={{42,24},{58,8}},
              lineColor={0,0,255},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid),
            Ellipse(
              extent={{-54,-8},{-38,-24}},
              lineColor={0,0,255},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid),
            Ellipse(
              extent={{12,-8},{28,-24}},
              lineColor={0,0,255},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid),
            Ellipse(
              extent={{-22,-40},{-6,-56}},
              lineColor={0,0,255},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid),
            Ellipse(
              extent={{12,-40},{28,-56}},
              lineColor={0,0,255},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid),
            Ellipse(
              extent={{42,-40},{58,-56}},
              lineColor={0,0,255},
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid),
            Bitmap(extent={{-56,130},{60,-6}}, fileName=
                  "C:/Users/Administrator/Desktop/h.png")}),
        Diagram(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-160,-100},{160,100}},
            grid={2,2},
            initialScale=0.5), graphics),
        Documentation(info="<html>
<p>
Model of a single U-tube borehole heat exchanger. 
The borehole heat exchanger is vertically discretized into <i>n<sub>seg</sub></i>
elements of height <i>h=h<sub>Bor</sub>&frasl;n<sub>seg</sub></i>.
Each segment contains a model for the heat transfer in the borehole, 
for heat transfer in the soil and for the far-field boundary condition.
</p>
<p>
The heat transfer in the borehole is computed using a convective heat transfer coefficient
that depends on the fluid velocity, a heat resistance between the two pipes, and
a heat resistance between the pipes and the circumference of the borehole.
The heat capacity of the fluid, and the heat capacity of the grout, is taken into account.
All thermal mass is assumed to be at the two bulk temperatures of the down-flowing 
and up-flowing fluid.
</p>
<p>
The heat transfer in the soil is computed using transient heat conduction in cylindrical
coordinates for the spatial domain <i>r<sub>bor</sub> &le; r &le; r<sub>ext</sub></i>. 
In the radial direction, the spatial domain is discretized into 
<i>n<sub>hor</sub></i> segments with uniform material properties.
Thermal properties can be specified separately for each horizontal layer.
The vertical heat flow is assumed to be zero, and there is assumed to be 
no ground water flow. 
</p>
<p>
The far-field temperature, i.e., the temperature at the radius 
<i>r<sub>ext</sub></i>, is computed using a power-series solution
to a line-source heat transfer problem. This temperature boundary condition
is updated every <i>t<sub>sample</sub></i> seconds.
</p>
<p>
The initial far-field temperature <i>T<sub>ext,start</sub></i>, which
is the temperature of the soil at a radius <i>r<sub>ext</sub></i>,
is computed 
as a function of the depth <i>z &gt; 0</i>. 
For a depth between <i>0 &le; z &le; z<sub>0</sub></i>, the temperature
is set to <i>T<sub>ext,0,start</sub></i>. 
The value of <i>z<sub>0</sub></i> is a parameter with a default of 10 meters.
However, there is large variability in the depth where the undisturbed soil temperature
starts.
For a depth of <i>z<sub>0</sub> &le; z &le; h<sub>bor</sub></i>,
the temperature is computed as
</p>
<p align=\"center\" style=\"font-style:italic;\">
  T<sup>i</sup><sub>ext,start</sub> = T<sub>ext,0,start</sub> + (z<sup>i</sup> - z<sub>0</sub>)  dT &frasl; dz
</p>
<p>
with <i>i &isin; {1, ..., n<sub>ver</sub>}</i>,
where the temperature gradient <i>dT &frasl; dz &ge; 0</i> is a parameter.
As with <i>z<sub>0</sub></i>, there is large variability in 
<i>dT &frasl; dz &ge; 0</i>. The default value is set to <i>1</i> Kelvin per 100 meters.
For the temperature of the grout, the same equations are applied, with
<i>T<sub>ext,0,start</sub></i> replaced with
<i>T<sub>fil,0,start</sub></i>, and 
<i>T<sup>i</sup><sub>ext,start</sub></i> replaced with
<i>T<sup>i</sup><sub>fil,start</sub></i>. 
The default setting uses the same temperature for the soil and the filling material.
</p>
<h4>Implementation</h4>
<p>
Each horizontal layer is modeled using an instance of
<a href=\"modelica://Buildings.Fluid.HeatExchangers.Boreholes.BaseClasses.BoreholeSegment\">
Buildings.HeatExchangers.Fluid.Boreholes.BaseClasses.BoreholeSegment</a>.
This model is composed of the model
<a href=\"modelica://Buildings.Fluid.HeatExchangers.Boreholes.BaseClasses.HexInternalElement\">
Buildings.Fluid.HeatExchangers.Boreholes.BaseClasses.HexInternalElement</a> which computes
the heat transfer in the pipes and the borehole filling,
of the model
<a href=\"modelica://Buildings.HeatTransfer.Conduction.SingleLayerCylinder\">
Buildings.HeatTransfer.Conduction.SingleLayerCylinder</a> which computes
the heat transfer in the soil, and
of the model
<a href=\"modelica://Buildings.Fluid.HeatExchangers.Boreholes.BaseClasses.TemperatureBoundaryCondition\">
Buildings.Fluid.HeatExchangers.Boreholes.BaseClasses.TemperatureBoundaryCondition</a> which computes
the far-field temperature boundary condition.
</p>
</html>",     revisions="<html>
<ul>
<li>
October 8, 2013, by Michael Wetter:<br/>
Removed parameter <code>show_V_flow</code>.
</li>
<li>
September 27, 2013, by Michael Wetter:<br/>
Added missing <code>each</code> in propagation of material properties.
</li>
<li>
August 2011, by Pierre Vigouroux:<br/>
First implementation.
</li>
</ul>
</html>"));
    end UTube;

    package TestBore

      model UTube "Model that tests the borehole model"
        extends Modelica.Icons.Example;
       package Medium = Buildings.Media.ConstantPropertyLiquidWater;
        BoreholeField.UTube borHol(
          redeclare each package Medium = Medium,
          dp_nominal=10000,
          dT_dz=0.0015,
          m_flow_nominal=0.3,
          redeclare each parameter
            Buildings.HeatTransfer.Data.BoreholeFillings.Bentonite matFil,
          samplePeriod=604800/5,
          rExt=0.3,
          positionsBor=DataFiles.readCSVmatrix("positions.csv"),
          redeclare parameter Buildings.HeatTransfer.Data.Soil.Sandstone matSoi(
              d=1300),
          hBor=100,
          TExt0_start=283.15,
          TFil0_start=283.15) "Borehole heat exchanger" annotation (Placement(
              transformation(extent={{-16,-36},{16,-4}}, rotation=0)));
            inner Modelica.Fluid.System system
          annotation (Placement(transformation(extent={{-100,-100},{-80,-80}})));
        Buildings.Fluid.Sources.Boundary_ph sin(          redeclare package
            Medium =
              Medium, nPorts=1) "Sink"
          annotation (Placement(transformation(extent={{56,-30},{36,-10}})));
        Buildings.Fluid.Sources.MassFlowSource_T sou(
          redeclare package Medium = Medium,
          use_m_flow_in=true,
          T=298.15,
          nPorts=1) "Source"
          annotation (Placement(transformation(extent={{-50,-30},{-30,-10}})));
        Modelica.Blocks.Sources.Pulse pulse(
          period=365*86400,
          startTime=365*86400/4,
          amplitude=0.6)
          annotation (Placement(transformation(extent={{-90,-30},{-70,-10}})));
      equation
        connect(pulse.y, sou.m_flow_in)       annotation (Line(
            points={{-69,-20},{-60,-20},{-60,-12},{-50,-12}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(sou.ports[1], borHol.port_a) annotation (Line(
            points={{-30,-20},{-10,-20}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(borHol.port_b, sin.ports[1]) annotation (Line(
            points={{10,-20},{36,-20}},
            color={0,127,255},
            smooth=Smooth.None));
       annotation(__Dymola_Commands(file="modelica://Buildings/Resources/Scripts/Dymola/Fluid/HeatExchangers/Boreholes/Examples/UTube.mos"
              "Simulate and plot"),
                Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,
                  -100},{100,100}}), graphics),
                        Documentation(info="<html>
<p>
This example models a borehole heat exchanger with two pipes that are
symmetrically spaced. 
The simulation period is 5 years.
From the 4th to the 10th months, the mass flow source switches on the 
flow rate through the borehole. The leaving
water of the mass flow source is <i>25</i>&deg;C,
and the water that returns from the borehole is between 
<i>20.5</i>&deg;C
and 
<i>21.5</i>&deg;C.
</p>
</html>",       revisions="<html>
<ul>
<li>
September 27, 2013, by Michael Wetter:<br/>
Corrected <code>StopTime</code> annotation.
</li>
<li>
August 2011, by Pierre Vigouroux:<br/>
First implementation.
</li>
</ul>
</html>"),experiment(
            StopTime=157680000,
            Tolerance=1e-05));
      end UTube;

      model heatRespondTest "Model that tests the borehole model"
        extends Modelica.Icons.Example;
       package Medium = Buildings.Media.ConstantPropertyLiquidWater(T_default=Modelica.SIunits.Conversions.from_degC(10));
        BoreholeField.UTube borHol(
          redeclare each package Medium = Medium,
          dp_nominal=10000,
          dT_dz=0.0015,
          m_flow_nominal=0.3,
          redeclare each parameter
            Buildings.HeatTransfer.Data.BoreholeFillings.Bentonite matFil,
          hBor=50,
          samplePeriod=604800/5,
          positionsBor=DataFiles.readCSVmatrix("positions.csv"),
          redeclare parameter Buildings.HeatTransfer.Data.Soil.Sandstone matSoi(
              d=1300),
          rExt=0.8,
          TExt0_start=283.15,
          TFil0_start=283.15) "Borehole heat exchanger" annotation (Placement(
              transformation(extent={{-16,-36},{16,-4}}, rotation=0)));
            inner Modelica.Fluid.System system
          annotation (Placement(transformation(extent={{-100,-100},{-80,-80}})));
        Buildings.Fluid.Sources.Boundary_ph sin(          redeclare package
            Medium =
              Medium, nPorts=1) "Sink"
          annotation (Placement(transformation(extent={{116,-32},{96,-12}})));
        Modelica.Fluid.Machines.ControlledPump pump(use_m_flow_set=true, redeclare
            package Medium = Medium,
          m_flow_nominal=0.15,
          p_a_nominal=150000,
          p_b_nominal=100000)
          annotation (Placement(transformation(extent={{6,40},{-12,58}})));
        Modelica.Blocks.Sources.Constant const(k=0.15)
          annotation (Placement(transformation(extent={{-50,72},{-30,92}})));
        Buildings.Fluid.MixingVolumes.MixingVolume vol(
          nPorts=2,
          redeclare package Medium = Medium,
          V=0.2,
          m_flow_nominal=0.15)
                 annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=-90,
              origin={50,44})));
        Buildings.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
          annotation (Placement(transformation(extent={{70,68},{52,86}})));
        Modelica.Blocks.Sources.Constant const1(k=1000)
          annotation (Placement(transformation(extent={{114,68},{94,88}})));
      equation
        connect(borHol.port_b, sin.ports[1]) annotation (Line(
            points={{10,-20},{54,-20},{54,-22},{96,-22}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(const.y, pump.m_flow_set) annotation (Line(
            points={{-29,82},{1.5,82},{1.5,56.38}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(borHol.port_b, vol.ports[1]) annotation (Line(
            points={{10,-20},{40,-20},{40,46}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(pump.port_a, vol.ports[2]) annotation (Line(
            points={{6,49},{14,49},{14,48},{40,48},{40,42}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(prescribedHeatFlow.port, vol.heatPort) annotation (Line(
            points={{52,77},{48,77},{48,54},{50,54}},
            color={191,0,0},
            smooth=Smooth.None));
        connect(prescribedHeatFlow.Q_flow, const1.y) annotation (Line(
            points={{70,77},{88,77},{88,78},{93,78}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(pump.port_b, borHol.port_a) annotation (Line(
            points={{-12,49},{-30,49},{-30,46},{-74,46},{-74,-20},{-10,-20}},
            color={0,127,255},
            smooth=Smooth.None));
       annotation(__Dymola_Commands(file="modelica://Buildings/Resources/Scripts/Dymola/Fluid/HeatExchangers/Boreholes/Examples/UTube.mos"
              "Simulate and plot"),
                Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},
                  {120,100}}),       graphics),
                        Documentation(info="<html>
<p>
This example models a borehole heat exchanger with two pipes that are
symmetrically spaced. 
The simulation period is 5 years.
From the 4th to the 10th months, the mass flow source switches on the 
flow rate through the borehole. The leaving
water of the mass flow source is <i>25</i>&deg;C,
and the water that returns from the borehole is between 
<i>20.5</i>&deg;C
and 
<i>21.5</i>&deg;C.
</p>
</html>",       revisions="<html>
<ul>
<li>
September 27, 2013, by Michael Wetter:<br/>
Corrected <code>StopTime</code> annotation.
</li>
<li>
August 2011, by Pierre Vigouroux:<br/>
First implementation.
</li>
</ul>
</html>"),experiment(
            StopTime=259200,
            Tolerance=1e-5),
          Icon(coordinateSystem(extent={{-100,-100},{120,100}})));
      end heatRespondTest;
    end TestBore;

    package positionsData
      function squareData
        input Integer   nPos = 2;
        input Real   space = 5;

       //   String fileName= "positions.csv";
       //   String headers[2]={"x","y"};

       //   constant String fileNames= fileName;

        output Real positions[nPos*nPos,2];
      protected
        Integer k;
      algorithm
        k:=1;
        for i in 1:nPos loop
        for j in 1:nPos loop
          positions[k,1]:= (i-1)*space;
          positions[k,2]:= (j-1)*space;
          k:=k+1;
        end for;
        end for;

        annotation ();
      end squareData;

      model squaresDataGFunctionWrite
         parameter String  fileName= "positions.csv";
         //   parameter String  headers[2]={"x","y"};
         parameter String  headers[2]={"x","y"};
          parameter Integer   nPos = 8;
        parameter Real   space = 5;

        // second parts
        parameter Real positionsBor[nPos*nPos,2]=BoreholeField.positionsData.squareData(nPos=nPos, space=space);
        parameter Real coeff[BoreholeField.Constants.N]=BoreholeField.coeffEi(
          r_b=rExt,
          H=hBor,
          positions=positionsBor);
          parameter Real coeffA[BoreholeField.Constants.N](each start=0);
         parameter Real coeff2[BoreholeField.Constants.N,2]=[coeff,coeffA];

        parameter Real rExt=0.3;
        parameter Real hBor=100;
           parameter String  fileNameG= "gFunction.csv";
         parameter String  headersG[1]={"gFunction"};

      //protected
       // Real pData[:,2];
      initial algorithm
      //  positionsBor:=BoreholeField.positionsData.squareData(nPos=nPos, space=space);

        DataFiles.writeCSVmatrix(
          fileName,
          headers,
          positionsBor,
          ",");
        DataFiles.writeCSVmatrix(
          fileNameG,
          headers,
          coeff2,
          ",");
      /*
  coeff:=BoreholeField.coeffEi(
    r_b=rExt,
    H=hBor,
    positions=positionsBor);*/
      end squaresDataGFunctionWrite;

      model squaresDataWrite
         parameter String  fileName= "positions.csv";
         //   parameter String  headers[2]={"x","y"};
         parameter String  headers[2]={"x","y"};
          parameter Integer   nPos = 1;
        parameter Real   space = 8;

      //protected
       // Real pData[:,2];
      initial algorithm

        DataFiles.writeCSVmatrix(
          fileName,
          headers,
          BoreholeField.positionsData.squareData(nPos=nPos, space=space),
          ",");
      end squaresDataWrite;

      model read
         parameter Real ac[:]=(DataFiles.readCSVmatrix("gFunction.csv"))[:,1];
      equation

      end read;
    end positionsData;
    annotation (
      conversion(noneFromVersion=""));
  end BoreholeField;

  package gshpSystem
    package BaseClasses
      partial model PartialElectric
        "Base class for test model of chiller electric EIR"
        import BoreholeField;
       package Medium1 = Buildings.Media.ConstantPropertyLiquidWater
          "Medium model";
       package Medium2 = Buildings.Media.ConstantPropertyLiquidWater
          "Medium model";

        parameter Modelica.SIunits.Power P_nominal=-per1.QEva_flow_nominal/per1.COP_nominal
          "Nominal compressor power (at y=1)";
        parameter Modelica.SIunits.TemperatureDifference dTEva_nominal=10
          "Temperature difference evaporator inlet-outlet";
        parameter Modelica.SIunits.TemperatureDifference dTCon_nominal=10
          "Temperature difference condenser outlet-inlet";
        parameter Real COPc_nominal = 3 "Chiller COP";
        parameter Modelica.SIunits.MassFlowRate mEva_flow_nominal=
           per1.mEva_flow_nominal "Nominal mass flow rate at evaporator";
        parameter Modelica.SIunits.MassFlowRate mCon_flow_nominal=
           per1.mCon_flow_nominal "Nominal mass flow rate at condenser";

        replaceable Buildings.Fluid.Chillers.BaseClasses.PartialElectric CHR1(show_T=
              true)
          constrainedby Buildings.Fluid.Chillers.BaseClasses.PartialElectric(
          redeclare package Medium1 = Medium1,
          redeclare package Medium2 = Medium2,
          energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
          dp1_nominal=6000,
          dp2_nominal=6000) "Chiller model"
          annotation (Placement(transformation(extent={{0,0},{20,20}})));
        inner Modelica.Fluid.System system
          annotation (Placement(transformation(extent={{284,-202},{304,-182}})));
        replaceable parameter Buildings.Fluid.Chillers.Data.BaseClasses.Chiller
                                                                                per1
        constrainedby Buildings.Fluid.Chillers.Data.BaseClasses.Chiller
          "Base class for performance data"
          annotation (Placement(transformation(extent={{110,44},{124,58}})));
        replaceable parameter Buildings.Fluid.Chillers.Data.BaseClasses.Chiller
                                                                                per2
        constrainedby Buildings.Fluid.Chillers.Data.BaseClasses.Chiller
          "Base class for performance data"
          annotation (Placement(transformation(extent={{110,10},{124,24}})));
        Buildings.Fluid.FixedResistances.FixedResistanceDpM pCon_out(
          redeclare package Medium = Medium1,
          m_flow_nominal=mCon_flow_nominal,
          dp_nominal=6000) "Flow resistance"
          annotation (Placement(transformation(extent={{64,70},{84,90}})));
        Modelica.Blocks.Sources.BooleanExpression booleanExpression(y=true)
          annotation (Placement(transformation(extent={{-60,38},{-40,58}})));
        Modelica.Blocks.Sources.RealExpression TSet(y=273.15 + 7)
          annotation (Placement(transformation(extent={{-60,22},{-40,42}})));

        Buildings.Fluid.Movers.FlowMachine_m_flow ERP(
          redeclare package Medium = Medium1,
          m_flow_nominal=mCon_flow_nominal,
          dp(start=214992),
          filteredSpeed=false,
          addPowerToMedium=false) "Condenser water pump" annotation (Placement(
              transformation(
              extent={{11,11},{-11,-11}},
              rotation=180,
              origin={-65,83})));
        Buildings.Fluid.FixedResistances.FixedResistanceDpM DZSC2(
          redeclare package Medium = Medium1,
          m_flow_nominal=mCon_flow_nominal,
          dp_nominal=6000) "Flow resistance"
          annotation (Placement(transformation(extent={{-118,74},{-98,94}})));
        Buildings.Fluid.FixedResistances.FixedResistanceDpM pCon_in(
          redeclare package Medium = Medium1,
          m_flow_nominal=mCon_flow_nominal,
          dp_nominal=6000) "Flow resistance"
          annotation (Placement(transformation(extent={{-40,74},{-20,94}})));
        Modelica.Blocks.Sources.RealExpression mCon_flow(y=0.5*mCon_flow_nominal)
          annotation (Placement(transformation(extent={{-92,92},{-72,112}})));
          /*
  Buildings.Fluid.Storage.ExpansionVessel DYBS2(redeclare package Medium =
        Medium1, V_start=1)
    annotation (Placement(transformation(extent={{-72,159},{-46,186}})));*/
        Buildings.Fluid.FixedResistances.FixedResistanceDpM EJSQ(
          redeclare package Medium = Medium1,
          m_flow_nominal=mCon_flow_nominal,
          dp_nominal=6000) "Flow resistance"
          annotation (Placement(transformation(extent={{8,136},{-22,158}})));
        Buildings.Fluid.FixedResistances.FixedResistanceDpM EFSQ(
          redeclare package Medium = Medium1,
          m_flow_nominal=mCon_flow_nominal,
          dp_nominal=6000) "Flow resistance"
          annotation (Placement(transformation(extent={{106,70},{126,90}})));
        Buildings.Fluid.FixedResistances.FixedResistanceDpM pEva_out1(
          redeclare package Medium = Medium2,
          m_flow_nominal=mCon_flow_nominal,
          dp_nominal=6000) "Flow resistance"
          annotation (Placement(transformation(extent={{50,-84},{70,-64}})));
        Buildings.Fluid.FixedResistances.FixedResistanceDpM DZSC1(
          redeclare package Medium = Medium2,
          m_flow_nominal=mCon_flow_nominal,
          dp_nominal=6000) "Flow resistance"
          annotation (Placement(transformation(extent={{-122,-88},{-102,-68}})));
        Buildings.Fluid.FixedResistances.FixedResistanceDpM pEva_in1(
          redeclare package Medium = Medium2,
          m_flow_nominal=mCon_flow_nominal,
          dp_nominal=6000) "Flow resistance"
          annotation (Placement(transformation(extent={{-60,-88},{-40,-68}})));
        Modelica.Blocks.Sources.RealExpression mCon_flow1(y=mEva_flow_nominal)
          annotation (Placement(transformation(extent={{-100,-68},{-80,-48}})));
          /*
  Buildings.Fluid.Storage.ExpansionVessel DYBS1(redeclare package Medium =
        Medium2, V_start=1)
        annotation (Placement(transformation(extent={{-110,-187},{-84,-160}})));
        */
        Buildings.Fluid.FixedResistances.FixedResistanceDpM KFSQ(
          redeclare package Medium = Medium2,
          m_flow_nominal=mCon_flow_nominal,
          dp_nominal=6000) "Flow resistance"
          annotation (Placement(transformation(extent={{26,-136},{-4,-114}})));
        Buildings.Fluid.FixedResistances.FixedResistanceDpM KJSQ(
          redeclare package Medium = Medium2,
          m_flow_nominal=mCon_flow_nominal,
          dp_nominal=6000) "Flow resistance"
          annotation (Placement(transformation(extent={{80,-84},{100,-64}})));
        Buildings.Fluid.FixedResistances.Pipe pip(redeclare package Medium = Medium2,
          thicknessIns=0.0001,
          lambdaIns=0.0001,
          m_flow_nominal=mEva_flow_nominal,
          nSeg=3,
          length=3)                                                                   annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=-90,
              origin={126,-106})));
        Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow
          prescribedHeatFlow
          annotation (Placement(transformation(extent={{170,-114},{152,-96}})));
        Modelica.Blocks.Sources.RealExpression heat(y=420000)
          annotation (Placement(transformation(extent={{144,-68},{164,-48}})));
        Buildings.Fluid.Actuators.Valves.TwoWayLinear V5(
          redeclare package Medium = Medium2,
          dpValve_nominal=20902,
          dpFixed_nominal=14930 + 89580,
          y_start=1,
          filteredOpening=false,
          m_flow_nominal=mEva_flow_nominal)
          "Control valve for chilled water leaving from chiller" annotation (
            Placement(transformation(
              extent={{-7,7},{7,-7}},
              rotation=0,
              origin={33,-73})));
        Buildings.Fluid.Actuators.Valves.TwoWayLinear V1(
          redeclare package Medium = Medium2,
          dpValve_nominal=20902,
          dpFixed_nominal=14930 + 89580,
          y_start=1,
          filteredOpening=false,
          m_flow_nominal=mEva_flow_nominal)
          "Control valve for chilled water leaving from chiller" annotation (
            Placement(transformation(
              extent={{-7,7},{7,-7}},
              rotation=0,
              origin={-23,-79})));
        Buildings.Fluid.Actuators.Valves.TwoWayLinear V3(
          redeclare package Medium = Medium1,
          dpValve_nominal=20902,
          dpFixed_nominal=14930 + 89580,
          y_start=1,
          filteredOpening=false,
          m_flow_nominal=mEva_flow_nominal)
          "Control valve for chilled water leaving from chiller" annotation (
            Placement(transformation(
              extent={{-7,7},{7,-7}},
              rotation=0,
              origin={-5,83})));
        Buildings.Fluid.Actuators.Valves.TwoWayLinear V7(
          redeclare package Medium = Medium1,
          dpValve_nominal=20902,
          dpFixed_nominal=14930 + 89580,
          y_start=1,
          filteredOpening=false,
          m_flow_nominal=mEva_flow_nominal)
          "Control valve for chilled water leaving from chiller" annotation (
            Placement(transformation(
              extent={{-7,7},{7,-7}},
              rotation=0,
              origin={45,81})));
        Buildings.Fluid.Actuators.Valves.TwoWayLinear V2(
          redeclare package Medium = Medium2,
          dpValve_nominal=20902,
          dpFixed_nominal=14930 + 89580,
          y_start=1,
          filteredOpening=false,
          m_flow_nominal=mEva_flow_nominal)
          "Control valve for chilled water leaving from chiller" annotation (
            Placement(transformation(
              extent={{-7,7},{7,-7}},
              rotation=0,
              origin={-23,-35})));
        Buildings.Fluid.Actuators.Valves.TwoWayLinear V6(
          redeclare package Medium = Medium2,
          dpValve_nominal=20902,
          dpFixed_nominal=14930 + 89580,
          y_start=1,
          filteredOpening=false,
          m_flow_nominal=mEva_flow_nominal)
          "Control valve for chilled water leaving from chiller" annotation (
            Placement(transformation(
              extent={{-7,7},{7,-7}},
              rotation=0,
              origin={35,-37})));
        Modelica.Blocks.Math.BooleanToReal scon "Contorl signal for chiller"
          annotation (Placement(transformation(extent={{-204,4},{-184,24}})));
        Modelica.Blocks.MathBoolean.Not nor1
          annotation (Placement(transformation(extent={{-218,-40},{-210,-32}})));
        Modelica.Blocks.Sources.BooleanExpression swh(y=true)
          annotation (Placement(transformation(extent={{-258,-30},{-238,-10}})));
        Modelica.Blocks.Math.BooleanToReal wcon "Contorl signal for chiller"
          annotation (Placement(transformation(extent={{-204,-46},{-184,-26}})));
        Buildings.Fluid.Actuators.Valves.TwoWayLinear V4(
          redeclare package Medium = Medium1,
          dpValve_nominal=20902,
          dpFixed_nominal=14930 + 89580,
          y_start=1,
          filteredOpening=false,
          m_flow_nominal=mEva_flow_nominal)
          "Control valve for chilled water leaving from chiller" annotation (
            Placement(transformation(
              extent={{-7,7},{7,-7}},
              rotation=0,
              origin={-5,117})));
        Buildings.Fluid.Actuators.Valves.TwoWayLinear V8(
          redeclare package Medium = Medium1,
          dpValve_nominal=20902,
          dpFixed_nominal=14930 + 89580,
          y_start=1,
          filteredOpening=false,
          m_flow_nominal=mEva_flow_nominal)
          "Control valve for chilled water leaving from chiller" annotation (
            Placement(transformation(
              extent={{-7,7},{7,-7}},
              rotation=0,
              origin={45,117})));
        Buildings.Fluid.Sources.FixedBoundary preSou2(redeclare package Medium
            = Medium1, nPorts=1)
          "Source for pressure and to account for thermal expansion of water"
          annotation (Placement(transformation(extent={{-188,114},{-166,136}})));
        Buildings.Fluid.Sources.FixedBoundary preSou1(
                                                     redeclare package Medium
            =                                                                   Medium2,
            nPorts=1)
          "Source for pressure and to account for thermal expansion of water"
          annotation (Placement(transformation(extent={{-158,-168},{-136,-146}})));

          EnergyStation.BoreholeField.UTube borHol(
          redeclare each package Medium = Medium1,
          dp_nominal=10000,
          dT_dz=0.0015,
          m_flow_nominal=0.3,
          redeclare each parameter
            Buildings.HeatTransfer.Data.BoreholeFillings.Bentonite matFil,
          redeclare parameter Buildings.HeatTransfer.Data.Soil.Sandstone matSoi,
          samplePeriod=604800/5,
          rExt=0.3,
          positionsBor=DataFiles.readCSVmatrix("positions.csv"),
          hBor=300,
          TExt0_start=283.15,
          TFil0_start=283.15) "Borehole heat exchanger" annotation (Placement(
              transformation(extent={{158,62},{190,94}}, rotation=0)));
      equation

        connect(booleanExpression.y,CHR1. on) annotation (Line(
            points={{-39,48},{-18,48},{-18,13},{-2,13}},
            color={255,0,255},
            smooth=Smooth.None));
        connect(CHR1.TSet, TSet.y) annotation (Line(
            points={{-2,7},{-28,7},{-28,32},{-39,32}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(DZSC2.port_b, ERP.port_a) annotation (Line(
            points={{-98,84},{-86,84},{-86,83},{-76,83}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(ERP.port_b, pCon_in.port_a) annotation (Line(
            points={{-54,83},{-46,83},{-46,84},{-40,84}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(EJSQ.port_b, DZSC2.port_a) annotation (Line(
            points={{-22,147},{-144,147},{-144,80},{-118,80},{-118,84}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(pCon_out.port_b, EFSQ.port_a) annotation (Line(
            points={{84,80},{106,80}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(KFSQ.port_b, DZSC1.port_a) annotation (Line(
            points={{-4,-125},{-130,-125},{-130,-78},{-122,-78}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(pEva_out1.port_b, KJSQ.port_a) annotation (Line(
            points={{70,-74},{80,-74}},
            color={0,127,255},
            smooth=Smooth.None));

        connect(KJSQ.port_b, pip.port_a) annotation (Line(
            points={{100,-74},{126,-74},{126,-96}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(pip.port_b, KFSQ.port_a) annotation (Line(
            points={{126,-116},{124,-116},{124,-136},{26,-136},{26,-125}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(pip.heatPort, prescribedHeatFlow.port) annotation (Line(
            points={{131,-106},{142,-106},{142,-105},{152,-105}},
            color={191,0,0},
            smooth=Smooth.None));
        connect(V5.port_b, pEva_out1.port_a) annotation (Line(
            points={{40,-73},{44,-73},{44,-74},{50,-74}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(pEva_in1.port_b, V1.port_a) annotation (Line(
            points={{-40,-78},{-36,-78},{-36,-80},{-30,-80},{-30,-79}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(V1.port_b, CHR1.port_a2) annotation (Line(
            points={{-16,-79},{-6,-80},{-6,-10},{20,-10},{20,4}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(pCon_in.port_b, V3.port_a) annotation (Line(
            points={{-20,84},{-16,84},{-16,83},{-12,83}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(V3.port_b, CHR1.port_a1) annotation (Line(
            points={{2,83},{4,83},{4,16},{0,16}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(pCon_out.port_a, V7.port_b) annotation (Line(
            points={{64,80},{58,80},{58,81},{52,81}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(V7.port_a, CHR1.port_b1) annotation (Line(
            points={{38,81},{34,81},{34,82},{20,82},{20,16}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(pEva_in1.port_b, V2.port_a) annotation (Line(
            points={{-40,-78},{-38,-78},{-38,-35},{-30,-35}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(V2.port_b, CHR1.port_a1) annotation (Line(
            points={{-16,-35},{-12,-35},{-12,16},{0,16}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(V6.port_b, pEva_out1.port_a) annotation (Line(
            points={{42,-37},{44,-37},{44,-74},{50,-74}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(CHR1.port_b1, V6.port_a) annotation (Line(
            points={{20,16},{28,16},{28,-37}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(CHR1.port_b2, V5.port_a) annotation (Line(
            points={{0,4},{0,-32},{26,-32},{26,-73}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(swh.y, scon.u) annotation (Line(
            points={{-237,-20},{-226,-20},{-226,14},{-206,14}},
            color={255,0,255},
            smooth=Smooth.None));
        connect(swh.y, nor1.u) annotation (Line(
            points={{-237,-20},{-226,-20},{-226,-36},{-219.6,-36}},
            color={255,0,255},
            smooth=Smooth.None));
        connect(nor1.y, wcon.u) annotation (Line(
            points={{-209.2,-36},{-206,-36}},
            color={255,0,255},
            smooth=Smooth.None));
        connect(scon.y, V1.y) annotation (Line(
            points={{-183,14},{-124,14},{-124,-106},{-23,-106},{-23,-87.4}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(V5.y, V1.y) annotation (Line(
            points={{33,-81.4},{33,-106},{-23,-106},{-23,-87.4}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(wcon.y, V2.y) annotation (Line(
            points={{-183,-36},{-90,-36},{-90,-48},{-23,-48},{-23,-43.4}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(V6.y, V2.y) annotation (Line(
            points={{35,-45.4},{35,-48},{-23,-48},{-23,-43.4}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(V4.port_b, CHR1.port_a2) annotation (Line(
            points={{2,117},{6,117},{6,28},{36,28},{36,4},{20,4}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(V8.port_a, CHR1.port_b2) annotation (Line(
            points={{38,117},{30,117},{30,120},{18,120},{18,32},{-8,32},{-8,4},{0,
                4}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(pCon_in.port_b, V4.port_a) annotation (Line(
            points={{-20,84},{-22,84},{-22,117},{-12,117}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(V8.port_b, pCon_out.port_a) annotation (Line(
            points={{52,117},{56,117},{56,112},{64,112},{64,80}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(scon.y, V3.y) annotation (Line(
            points={{-183,14},{-122,14},{-122,64},{-5,64},{-5,74.6}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(V7.y, V3.y) annotation (Line(
            points={{45,72.6},{45,64},{-5,64},{-5,74.6}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(wcon.y, V4.y) annotation (Line(
            points={{-183,-36},{-98,-36},{-98,108.6},{-5,108.6}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(V4.y, V8.y) annotation (Line(
            points={{-5,108.6},{45,108.6}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(preSou2.ports[1], DZSC2.port_a) annotation (Line(
            points={{-166,125},{-142,125},{-142,122},{-144,126},{-144,80},{-118,80},{-118,
                84}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(preSou1.ports[1], DZSC1.port_a) annotation (Line(
            points={{-136,-157},{-132,-157},{-132,-160},{-128,-160},{-128,-126},{-130,
                -125},{-130,-78},{-122,-78}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(EFSQ.port_b, borHol.port_a) annotation (Line(
            points={{126,80},{142,80},{142,78},{164,78}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(borHol.port_b, EJSQ.port_a) annotation (Line(
            points={{184,78},{222,78},{222,147},{8,147}},
            color={0,127,255},
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-360,
                  -220},{320,280}}),
                            graphics), Icon(coordinateSystem(extent={{-360,-220},{320,
                  280}})));
      end PartialElectric;
    end BaseClasses;




    model ElectricEIRT2dddd
      extends gshpSystem.BaseClasses.PartialElectric(
        redeclare Buildings.Fluid.Chillers.ElectricEIR CHR1(
          per=per1,
          T1_start=273.15 + 35,
          T2_start=273.15 + 10,
          show_T=true,
          m2_flow(fixed=false, start=20)),
        redeclare parameter
          Buildings.Fluid.Chillers.Data.ElectricEIR.ElectricEIRChiller_McQuay_WSC_471kW_5_89COP_Vanes
          per1 "CHR1, PSRHH-1201D , COOLC 480.4KW, HEATC 490.9KW",
        redeclare parameter
          Buildings.Fluid.Chillers.Data.ElectricEIR.ElectricEIRChiller_Trane_CGWD_207kW_3_99COP_None
          per2 "CHR2, PSRHH-Y0851R, COOLC 220KW, HEATC 234KW",
        borHol(dT_dz=0, hBor=100),
        TSet(y=273.15 + 3),
        mCon_flow(y=mCon_flow_nominal*0.35),
        pip(
          m_flow_nominal=mEva_flow_nominal,
          lambdaIns=0.01,
          T_start=300.15),
        heat(y=100000*sin(time/100000)),
        swh(y=mLoad.y > 100),
        mCon_flow1(y=2),
        system(allowFlowReversal=false),
        booleanExpression(y=if abs(heatP.y) > 1000 then true else false));
      Modelica.Blocks.Sources.RealExpression TSet4(y=5)
        annotation (Placement(transformation(extent={{-276,162},{-256,182}})));
      Modelica.Blocks.Continuous.LimPID PID1(
        controllerType=Modelica.Blocks.Types.SimpleController.PI,
        initType=Modelica.Blocks.Types.InitPID.InitialState,
        yMax=6,
        yMin=2,
        Ti=120,
        k=0.1)
        annotation (Placement(transformation(extent={{-230,188},{-210,208}})));
      Modelica.Blocks.Sources.RealExpression TSet5(y=(if swh.y then abs(CHR1.sta_a1.T
             - CHR1.sta_b1.T) else abs(CHR1.sta_a2.T - CHR1.sta_b2.T)))
        annotation (Placement(transformation(extent={{-276,204},{-256,224}})));
      Modelica.Blocks.Tables.CombiTable1D combiTable1D(
        tableOnFile=true,
        tableName="tabLoad",
        fileName="load.txt")
        annotation (Placement(transformation(extent={{162,-34},{182,-14}})));
      Modelica.Blocks.Sources.RealExpression TSet1(y=mod(time, 31536000))
        annotation (Placement(transformation(extent={{98,-30},{118,-10}})));
      Modelica.Blocks.Math.Product heatP
        annotation (Placement(transformation(extent={{210,-50},{230,-30}})));
      Modelica.Blocks.Sources.RealExpression consM(y=-1000)
        annotation (Placement(transformation(extent={{184,-86},{204,-66}})));
      Modelica.Blocks.Sources.RealExpression mLoad(y=if abs(heatP.y) > 1000
             then heatP.y else 0)
        annotation (Placement(transformation(extent={{160,-146},{180,-126}})));
      Buildings.Fluid.Movers.FlowMachine_m_flow ERP1(
        redeclare package Medium = Medium1,
        m_flow_nominal=mCon_flow_nominal,
        dp(start=214992),
        filteredSpeed=false,
        addPowerToMedium=false) "Condenser water pump" annotation (Placement(
            transformation(
            extent={{11,11},{-11,-11}},
            rotation=180,
            origin={-79,-83})));
    equation

      connect(TSet5.y, PID1.u_s) annotation (Line(
          points={{-255,214},{-232,214},{-232,198}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(TSet4.y, PID1.u_m) annotation (Line(
          points={{-255,172},{-220,172},{-220,186}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(PID1.y, ERP.m_flow_in) annotation (Line(
          points={{-209,198},{-58,198},{-58,96.2},{-65.22,96.2}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(TSet1.y, combiTable1D.u[1]) annotation (Line(
          points={{119,-20},{160,-20},{160,-24}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(combiTable1D.y[1], heatP.u1) annotation (Line(
          points={{183,-24},{208,-24},{208,-34}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(consM.y, heatP.u2) annotation (Line(
          points={{205,-76},{208,-76},{208,-46}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(mLoad.y, prescribedHeatFlow.Q_flow) annotation (Line(
          points={{181,-136},{170,-136},{170,-105}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(PID1.y, ERP1.m_flow_in) annotation (Line(
          points={{-209,198},{-209,-69.8},{-79.22,-69.8}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(DZSC1.port_b, ERP1.port_a) annotation (Line(
          points={{-102,-78},{-98,-78},{-98,-92},{-90,-92},{-90,-83}},
          color={0,127,255},
          smooth=Smooth.None));
      connect(ERP1.port_b, pEva_in1.port_a) annotation (Line(
          points={{-68,-83},{-68,-94},{-60,-94},{-60,-78}},
          color={0,127,255},
          smooth=Smooth.None));
    annotation (
    experiment(StopTime=14400), Diagram(coordinateSystem(preserveAspectRatio=false,
              extent={{-360,-220},{320,280}}), graphics));
    end ElectricEIRT2dddd;
    annotation ();
    model timet
        Real sd;
    equation
        sd= mod(time,450);
    end timet;
  end gshpSystem;
  annotation (uses(
      Buildings(version="1.6"),
      Modelica(version="3.2.1"),
      DataFiles(version="1.0.1"),
      BoreholeField(version="1")));
end EnergyStation;
