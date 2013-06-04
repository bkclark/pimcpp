---
layout: base
title: Constructing-a-density-matrix
---

Our task now is to construct a pair action between two He-4 particles.
squarer is a code by David Ceperley designed to perform matrix squaring
based on the specification of one of a number of classical pair
interactions. We also note the existence of a similar code, squarer++,
by Ken Esler.

Here, we will use squarer to generate the squared pair density matrix
and then we will format the output for use in PIMC++. Because this
involves several steps and can introduce frustrating errors, we are
providing a script to handle this process. The script is launched by
invoking

    generateMatrix

It will prompt you for the filename, which should be something like
"He.dm". To see what files have been generated once the script
completes, you can issue a command to display the most recently modified
files in the directory:

    ls -ltrh

Here, the relevant files are described.

We start with a file that contains the classical pair potential between
two particles:

     UNITS K A 
     TYPE He4 6.059615 
     TYPE He4 6.059615 
     GRID 250 LINEAR .2 8.0 
     SQUARER 10. 4 3 3 30 14 
     POTTAIL   -0.21715E-06  3  -6.00 -0.10124E+05  -8.00 -0.27612E+05 -10.00 -0.10179E+06
     RANK 2  250 1
     GRID 1 LINEAR .2 8.0 
     LABEL 1 r
     BEGIN potential 0
      0.153629521E+07  0.121764428E+07  0.992336768E+06  0.826575151E+06  0.700295964E+06
      0.601088635E+06  0.521030739E+06  0.454933193E+06  0.399319317E+06  0.351808764E+06
      0.310735332E+06  0.274904828E+06  0.243439233E+06  0.215675185E+06  0.191097242E+06
      0.169293676E+06  0.149926999E+06  0.132714182E+06  0.117413267E+06  0.103814220E+06
      0.917325762E+05  0.810049605E+05  0.714858306E+05  0.630450462E+05  0.555659862E+05
      0.489440375E+05  0.430853401E+05  0.379057112E+05  0.333297033E+05  0.292897631E+05
      0.257254715E+05  0.225828533E+05  0.198137468E+05  0.173752271E+05  0.152290800E+05
      0.133413211E+05  0.116817584E+05  0.102235945E+05  0.894306647E+04  0.781911988E+04
      0.683311469E+04  0.596856071E+04  0.521087985E+04  0.454719298E+04  0.396612909E+04
      0.345765480E+04  0.301292214E+04  0.262413296E+04  0.228441807E+04  0.198772992E+04
      0.172874708E+04  0.150278953E+04  0.130510159E+04  0.111588524E+04  0.953004566E+03
      0.812858391E+03  0.692333170E+03  0.588737573E+03  0.499745711E+03  0.423347895E+03
      0.357807947E+03  0.301626232E+03  0.253507639E+03  0.212333867E+03  0.177139449E+03
      0.147091007E+03  0.121469312E+03  0.996537657E+02  0.811089927E+02  0.653732379E+02
      0.520483441E+02  0.407910851E+02  0.313056762E+02  0.233372998E+02  0.166665104E+02
      0.111043998E+02  0.648841771E+01  0.267876156E+01 -0.444743897E+00 -0.298534230E+01
     -0.503167472E+01 -0.665979175E+01 -0.793489144E+01 -0.891282046E+01 -0.964137044E+01
     -0.101613973E+02 -0.105077877E+02 -0.107102932E+02 -0.107942500E+02 -0.107812002E+02
     -0.106894275E+02 -0.105344194E+02 -0.103292653E+02 -0.100850001E+02 -0.981089957E+01
     -0.951473526E+01 -0.920299380E+01 -0.888106542E+01 -0.855340622E+01 -0.822367758E+01
     -0.789486591E+01 -0.756938547E+01 -0.724916649E+01 -0.693573078E+01 -0.663025641E+01
     -0.633363301E+01 -0.604650915E+01 -0.576933262E+01 -0.550238487E+01 -0.524581026E+01
     -0.499964098E+01 -0.476381810E+01 -0.453843229E+01 -0.432367400E+01 -0.411919280E+01
     -0.392460844E+01 -0.373953183E+01 -0.356357012E+01 -0.339633082E+01 -0.323742509E+01
     -0.308647046E+01 -0.294309281E+01 -0.280692807E+01 -0.267762335E+01 -0.255483787E+01
     -0.243824350E+01 -0.232752516E+01 -0.222238098E+01 -0.212252235E+01 -0.202767379E+01
     -0.193757278E+01 -0.185196947E+01 -0.177062637E+01 -0.169331793E+01 -0.161983019E+01
     -0.154996029E+01 -0.148351605E+01 -0.142031551E+01 -0.136018649E+01 -0.130296609E+01
     -0.124850031E+01 -0.119664355E+01 -0.114725821E+01 -0.110021428E+01 -0.105538895E+01
     -0.101266620E+01 -0.971936422E+00 -0.933096125E+00 -0.896047544E+00 -0.860698340E+00
     -0.826961295E+00 -0.794754022E+00 -0.763998686E+00 -0.734621753E+00 -0.706553735E+00
     -0.679728963E+00 -0.654085368E+00 -0.629564269E+00 -0.606110184E+00 -0.583670639E+00
     -0.562196002E+00 -0.541639313E+00 -0.521956134E+00 -0.503104403E+00 -0.485044298E+00
     -0.467738106E+00 -0.451150109E+00 -0.435246463E+00 -0.419995095E+00 -0.405365605E+00
     -0.391329166E+00 -0.377858438E+00 -0.364927487E+00 -0.352511703E+00 -0.340587725E+00
     -0.329133378E+00 -0.318127600E+00 -0.307550384E+00 -0.297382722E+00 -0.287606547E+00
     -0.278204683E+00 -0.269160798E+00 -0.260459357E+00 -0.252085578E+00 -0.244025397E+00
     -0.236265423E+00 -0.228792907E+00 -0.221595705E+00 -0.214662250E+00 -0.207981519E+00
     -0.201543005E+00 -0.195336691E+00 -0.189353025E+00 -0.183582896E+00 -0.178017611E+00
     -0.172648872E+00 -0.167468762E+00 -0.162469719E+00 -0.157644521E+00 -0.152986272E+00
     -0.148488380E+00 -0.144144548E+00 -0.139948753E+00 -0.135895240E+00 -0.131978504E+00
     -0.128193277E+00 -0.124534522E+00 -0.120997416E+00 -0.117577345E+00 -0.114269890E+00
     -0.111070819E+00 -0.107976080E+00 -0.104981789E+00 -0.102084227E+00 -0.992798276E-01
     -0.965651731E-01 -0.939369863E-01 -0.913921244E-01 -0.889275729E-01 -0.865404399E-01
     -0.842279502E-01 -0.819874405E-01 -0.798163544E-01 -0.777122372E-01 -0.756727317E-01
     -0.736955741E-01 -0.717785895E-01 -0.699196883E-01 -0.681168624E-01 -0.663681815E-01
     -0.646717903E-01 -0.630259045E-01 -0.614288085E-01 -0.598788520E-01 -0.583744473E-01
     -0.569140671E-01 -0.554962413E-01 -0.541195550E-01 -0.527826463E-01 -0.514842037E-01
     -0.502229644E-01 -0.489977122E-01 -0.478072756E-01 -0.466505256E-01 -0.455263747E-01
     -0.444337747E-01 -0.433717152E-01 -0.423392222E-01 -0.413353566E-01 -0.403592128E-01
     NDERIV  2

We will now look at this file in detail. The second and third line
specify the type of this density matrix (in this case He4) and its
respective lambda (\\(\lambda =   \frac{\hbar^2}{2m}\\)). The third line
specifies the grid we are putting our potential on. The first number
after "GRID" specifies the number of points on our grid. The word
"LINEAR" specifies that it will be a linear grid and the remaining two
numbers specify the range of the grid (from 0.2 to 8.0). The next line

     SQUARER 10. 4 3 3 30 14 

is the key line in telling squarer what to do. It specifies (lowest
temperature reached, number of tau kept in file, dimension, something
else, number of plane waves used in the construction, number of
squarings done to reach the lowest temperature).

We will look at each of these in turn. **Lowest Temperature Reached**:
This is (the inverse of) the lowest tau in the file and would typically
be the time step \\(\tau\\) that you intend to be using in the simulation
(or possibly a little lower so you can do a convergence study in time
step).

**Number of tau kept in file:** It is important to keep at least as many
tau as there are levels in your BisectionBlock move (the BisectionBlock
move uses these larger tau's to sample more effectively) (it is not
strictly illegal to use the BisectionBlock without these tau's but will
degrade your sampling effeciency and should not typically be done).

**Dimension**: Self explanatory

**Plane waves used**: Squarer decomposes the potential into plane waves
internally (and then converts back before writing the output). You
should choose a number large enough so that your result is converged.
Squarer will tell you at the end of the run how many plane waves have
been used. If it turns out that it used exactly the number that you
chose, then you probably need to choose a bigger number.

**Squarings done to reach lowest temperature:** The temperature that the
squarer program starts at is 2\^(squarings done). this means that you
want to ensure that your lowest temperature time 2\^(squarings done) is
well approximated by the primitive approximation.

The next line

      POTTAIL   -0.21715E-06  3  -6.00 -0.10124E+05  -8.00 -0.27612E+05 -10.00 -0.10179E+06 

is not used by squarer but is supposed to encode information about the
way in which you shifted your potential so that it wasn't discontinuous
at the half box (Note: if you haven't done this, there is a reasonable
chance you've made a mistake).

The final lines specify the potential on the grid. The data is placed in
five columns, but these five columns should just be treated as a flat
set of numbers on which the potential is specified.

Squarer takes the parameters and pair potential specified in the .dm
file and adds to the .dm file data about the squared density matrices.

For use in PIMC++, we need files with extension .Sampling.in and
.PairAction. The generateMatrix script will produce these.

The .PairAction file should look like this:

    string SamplingTableFile="He4.95.Sampling.in";

    Section (Fits)
    {
      string Type = "DavidFit";
      int NumOffDiagonalTerms=2;
      Section (Particle1)
      {
        string Name = "He";
        double lambda = 6.059615;
        int Ndim = 3;
      }

      Section (Particle2)
      {
        string Name = "He";
        double lambda = 6.059615;
        int Ndim = 3;
      }

      string Daviddmfile="He4.95.dm";
    }

where it has the appropriate Type names (that are used in the rest of
the input file), the appropriate path to the SamplingTableFile and the
appropriate dmFile.
