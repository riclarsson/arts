/* Copyright (C) 2000 Stefan Buehler <sbuehler@uni-bremen.de>,
                      Axel von Engeln <engeln@uni-bremen.de>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

/*!
  \file   species_data.cc
  \brief  Implementation of function define_species_data().

  This file contains the definition of this function and nothing
  else. You can add new records here if you want to extend the
  capability of ARTS.

  \author Stefan Buehler, Axel von Engeln
  \date 2000-08-10 */

#include "arts.h"
#include "make_array.h"
#include "absorption.h"


/*! The lookup information for all the different species. */
ARRAY<SpeciesRecord> species_data;


/*! \name Some #defines for better readability */
//@{ 
#define NAME(x) x				
#define DEGFR(x) x				
#define ISOTOPES     make_array<IsotopeRecord>	
#define REC	     IsotopeRecord		
#define TAGS	     make_array<int>		
//@} 


/*!
  \brief Define species lookup data.

  <h1>Source for entries</h1>
  <!------------------------>

  <dl>
  <dt> Molecule Name:
  <dd> Generally Hitran convention, unless otherwise indicated.

  <dt> Degrees of freedom: 
  <dd> Generally taken from the forward_4_96 release,
       file glob_def.c, unless otherwise indicated.

  <dt> Isotopes: 
  <dd> Generally Hitran convention, e.g., for the main H2O
       isotope: 161, thus the last digit of the sum of the
       neutrons and protons of the individual atoms is
       taken. Nevertheless, molecules like H2CO have only 3
       digits in Hitran, only one for the two H atoms. This leads
       to problems with the jpl catalogue, if this includes DHCO
       and DDCO. For these cases, a different scheme is taken as
       indicated.

  <dt> Abundance: 
  <dd> Generally Hitran convention, unless otherwise indicated.
       A number code in the header gives the source of the abundancies 
       for each isotope of each species:
       
       <table>
       <tr>
       <td> 1: <td> hitran 96 abundance, taken from cd: file
          	    software/generic/tables_96.txt (default, even if jpl and
          	    hitran abundancies are available).
       <tr>
       <td >2: <td> jpl abundance, taken from the documentation coming
          	    along with the catalogue. latest catalogue version
          	    extracted 27.07.00, can be found at
          	    /pool/lookup/jpl/cat7_00/doc/d<tag_nr>.cat 
       <tr>
       <td> 3: <td> jpl abundance is multiplied with the maximum abundance
          	    of this species found in hitran. is only performed when
          	    abundancies of 1 were found in the jpl catalogue.
       </table>

  <dt> Mass: 
  <dd> Rounded atomic weight. (SAB: From looking at the table in forward_4_96,
       glob_def.c, the difference between the actual mass and the mass
       simply estimated from the Atom number is only 0.001. This seems
       not worth the trouble. Anyway, the field for the mass is
       there. Should anybody feel like adding the true numbers, just go
       ahead.)

  <dt >MY-tags:
  <dd> Extracted from file glob_def.c.

  <dt> HI-tags:
  <dd> From hitran 96, software/generic/tables_96.txt.

  <dt> JPL-tags:
  <dd> Collected in file 
       /pool/lookup/jpl/cat7_00/abundancies/tag_species.jpl, taken from
       last issue of jpl catalogue (7/00).
  </dl>

  Some more information can be found at /pool/lookup/jpl/cat7_00/abundancies,
  where the idl script that reads/converts the abundances is located.

  08.08.00 AvE

  \author Stefan Buehler, Axel von Engeln  

  <b>History</b>
  <dl>
  <dt> 2000-08-08 Axel von_Engeln:
  <dd> Added all JPL/HITRAN species and documentation.
  </dl>
*/
void define_species_data()
{
  extern ARRAY<SpeciesRecord> species_data;

  // Initialize to zero, just in case:
  species_data.clear();

  /* Here's an empty template record entry:

  species_data.push_back
    ( SpeciesRecord
      ( NAME("H2O"),
	DEGFR(3),
	ISOTOPES
	(//   Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC( ""	,		,	,	,	,TAGS() ),
	 REC( ""	,		,	,	,	,TAGS() )
	 ) ) );

  For good looks, keep the commas on the marks!

  */


  // H2O
  // Abundancies: 1 1 1 1 3 3
  species_data.push_back
    ( SpeciesRecord
      ( NAME("H2O"),
	DEGFR(3),
	ISOTOPES
	(//   Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC( "161"	,0.99731702	,18.	,11	,11	,TAGS(18003, 18005) ),
	 REC( "181"	,0.00199983	,20.	,12	,12	,TAGS(20003) ),
	 REC( "171"	,0.00037200	,19.	,13	,13	,TAGS(19003) ),
	 REC( "162"	,0.00031069	,19.	,14	,14	,TAGS(19002) ),
	 REC( "182"	,6.1070746E-07	,21.	,-1	,-1	,TAGS(21001) ),
	 REC( "262"	,2.2430204E-08	,20.	,-1	,-1	,TAGS(20001) )
	 ) ) );

  // CO2 
  // (missing mainly in JPL, latest version (7/00) includes some isotopes)
  // Degrees of freedom from Schanda:`Physical Fundamentals of Remote Sensing'
  // Abundancies: 1 1 1 1 1 1 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("CO2"),
	DEGFR(2),
	ISOTOPES
	(//   Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC( "626"	,0.98420	,44.	,21	,21	,TAGS() ),
	 REC( "636"	,0.01106	,45.	,22	,22	,TAGS() ),
	 REC( "628"	,0.0039471	,46.	,23	,23	,TAGS(46013) ),
	 REC( "627"	,0.000734	,45.	,24	,24	,TAGS(45012) ),
	 REC( "638"	,0.00004434	,47.	,25	,25	,TAGS() ),
	 REC( "637"	,0.00000825	,46.	,26	,26	,TAGS() ),
	 REC( "828"	,0.0000039573	,48.	,27	,27	,TAGS() ),
	 REC( "728"	,0.00000147	,47.	,28	,28	,TAGS() )
	 ) ) );
  

  // O3
  // Abundancies: 1 1 1 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("O3"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("666"	,0.992901	,48.	,31	,31	,TAGS(48004, 48005, 48006, 48007, 48008)),
	 REC("668"	,0.00398194	,50.	,32	,32	,TAGS(50004, 50006)),
	 REC("686"	,0.00199097	,50.	,33	,33	,TAGS(50003, 50005)),
	 REC("667"	,0.000740	,49.	,34	,34	,TAGS(49002)),
	 REC("676"	,0.000370	,49.	,35	,35	,TAGS(49001))
	 ) ) );


  // N2O
  // Abundancies: 1 1 1 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("N2O"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("446"	,0.990333	,44.	,41	,41	,TAGS(44004, 44009, 44012)),
	 REC("456"	,0.0036409	,45.	,42	,42	,TAGS(45007)),
	 REC("546"	,0.0036409	,45.	,43	,43	,TAGS(45008)),
	 REC("448"	,0.00198582	,46.	,44	,44	,TAGS(46007)),
	 REC("447"	,0.000369	,45.	,-1	,45	,TAGS() )
	 ) ) );

  // CO
  // Abundancies: 1 1 1 1 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("CO"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("26"	,0.98654	,28.	,51	,51	,TAGS(28001)),
	 REC("36"	,0.01108	,29.	,52	,52	,TAGS(29001)),
	 REC("28"	,0.0019782	,30.	,53	,53	,TAGS(30001)),
	 REC("27"	,0.000368	,29.	,-1	,54	,TAGS(29006)),
	 REC("38"	,0.00002222	,31.	,-1	,55	,TAGS() ),
	 REC("37"	,0.00000413	,30.	,-1	,56	,TAGS() )
	 ) ) );

  // CH4
  // Degrees of freedom: jpl catalogue
  // Abundancies: 1 1 1
  // Note: - jpl abundancies for tag 17003: 0.00014996848
  //       - CH4 is in official mytran list (6), but does not 
  //         seem to be included for calculation, as given
  //         by table tag_table in file glob_def.c
  species_data.push_back
    ( SpeciesRecord
      ( NAME("CH4"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("211"	,0.98827	,16.	,-1	,61	,TAGS()),
	 REC("311"	,0.01110	,17.	,-1	,62	,TAGS()),
	 REC("212"	,0.00061575	,17.	,-1	,63	,TAGS(17003))
	 ) ) );

  // O2
  // Abundancies: 1 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("O2"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("66"	,0.995262	,32.	,71	,71	,TAGS(32001, 32002, 32005)),
	 REC("68"	,0.00399141	,34.	,72	,72	,TAGS(34001)),
	 REC("67"	,0.000742	,33.	,73	,73	,TAGS(33002))
	 ) ) );

  // NO
  // Abundancies: 1 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("NO"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("46"	,0.993974	,30.	,81	,81	,TAGS(30008)),
	 REC("56"	,0.0036543	,31.	,-1	,82	,TAGS() ),
	 REC("48"	,0.00199312	,32.	,-1	,83	,TAGS() )
	 ) ) );

  // SO2
  // Abundancies: 1 1 2 2
  species_data.push_back
    ( SpeciesRecord
      ( NAME("SO2"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("626"	,0.94568	,64.	,91	,91	,TAGS(64002, 64005)),
	 REC("646"	,0.04195	,66.	,-1	,92	,TAGS(66002)),
	 REC("636"	,0.0074989421	,65.	,-1	,-1	,TAGS(65001)),
	 REC("628"	,0.0020417379	,66.	,-1	,-1	,TAGS(66004))
	 ) ) );

  // NO2
  // Abundancies: 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("NO2"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("646"	,0.991616	,46.	,101	,101	,TAGS(46006))
	 ) ) );

  // NH3
  // Abundancies: 1 1 3
  species_data.push_back
    ( SpeciesRecord
      ( NAME("NH3"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("4111"	,0.9958715	,17.	,111	,111	,TAGS(17002, 17004)),
	 REC("5111"	,0.0036613	,18.	,112	,112	,TAGS(18002)),
	 REC("4112"	,0.00044792294	,18.	,-1	,-1	,TAGS(18004))
	 ) ) );

  // HNO3
  // Abundancies: 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("HNO3"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("146"	,0.989110	,63.	,121	,121	,TAGS(63001, 63002, 63003, 63004, 63005, 63006))
	 ) ) );

  // OH
  // Abundancies: 1 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("OH"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("61"	,0.997473	,17.	,131	,131	,TAGS(17001)),
	 REC("81"	,0.00200014	,19.	,132	,132	,TAGS(19001)),
	 REC("62"	,0.00015537	,18.	,133	,133	,TAGS(18001))
	 ) ) );

  // HF
  // Abundancies: 1 3
  species_data.push_back
    ( SpeciesRecord
      ( NAME("HF"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("19"	,0.99984425	,20.	,141	,141	,TAGS(20002)),
	 REC("29"	,0.00014994513	,21.	,-1	,-1	,TAGS(21002))
	 ) ) );

  // HCl
  // Abundancies: 1 1 2 2
  species_data.push_back
    ( SpeciesRecord
      ( NAME("HCl"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("15"	,0.757587	,36.	,151	,151	,TAGS(36001)),
	 REC("17"	,0.242257	,38.	,152	,152	,TAGS(38001)),
	 REC("25"	,0.00011324004	,37.	,-1	,-1	,TAGS(37001)),
	 REC("27"	,3.6728230E-05	,39.	,-1	,-1	,TAGS(39004))
	 ) ) );

  // HBr
  // Abundancies: 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("HBr"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("19"	,0.50678	,80.	,161	,161	,TAGS(80001)),
	 REC("11"	,0.49306	,82.	,162	,162	,TAGS(82001))
	 ) ) );

  // HI
  // Degrees of freedom: guessed, since it seems to be linear
  // Abundancies: 1
  // Note: HI is in official mytran list (17), but does not 
  //       seem to be included for calculation, as given
  //       by table tag_table in file glob_def.c

  species_data.push_back
    ( SpeciesRecord
      ( NAME("HI"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("17"	,0.99984425	,128.	,-1	,171	,TAGS( ))
	 ) ) );

  // ClO
  // Abundancies: 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("ClO"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("56"	,0.75591	,51.	,181	,181	,TAGS(51002, 51003)),
	 REC("76"	,0.24172	,53.	,182	,182	,TAGS(53002, 53006))
	 ) ) );

  // OCS
  // Abundancies: 1 1 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("OCS"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("622"	,0.93739	,60.	,191	,191	,TAGS(60001)),
	 REC("624"	,0.04158	,62.	,192	,192	,TAGS(62001)),
	 REC("632"	,0.01053	,61.	,193	,193	,TAGS(61001)),
	 REC("822"	,0.0018797	,62.	,194	,194	,TAGS(62002))
	 ) ) );

  // H2CO
  // Abundancies: 1 1 1 3 3
  // Note: the isotope names differ from hitran convention, since the jpl catalogue has 
  //       isotopes HHCO, HDCO, DDCO.
  //       hitran convention  --  new convention  -- jpl species
  //              126                 1126             H2CO
  //              136                 1136             H2C-13-O
  //              128                 1128             H2CO-18
  //               -                  1226             HDCO
  //               -                  2226             D2CO
  species_data.push_back
    ( SpeciesRecord
      ( NAME("H2CO"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("1126"	,0.98624	,30.	,201	,201	,TAGS(30004)),
	 REC("1136"	,0.01108	,31.	,202	,202	,TAGS(31002)),
	 REC("1128"	,0.0019776	,32.	,203	,203	,TAGS(32004)),
	 REC("1226"	,0.00029578940	,31.	,-1	,-1	,TAGS(31003)),
	 REC("2226"	,2.2181076E-08	,32.	,-1	,-1	,TAGS(32006))
	 ) ) );

  // HOCl
  // Abundancies: 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("HOCl"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("165"	,0.75579	,52.	,211	,211	,TAGS(52006)),
	 REC("167"	,0.24168	,54.	,212	,212	,TAGS(54005))
	 ) ) );

  // N2
  // Degrees of freedom: guessed, since it seems to be linear
  // Abundancies: 1
  // Note: N2 is in official mytran list (22), but does not 
  //       seem to be included for calculation, as given
  //       by table tag_table in file glob_def.c
  species_data.push_back
    ( SpeciesRecord
      ( NAME("N2"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("44"	,0.9926874	,28.	,-1	,221	,TAGS( ))
	 ) ) );

  // HCN
  // Abundancies: 1 1 1 3
  species_data.push_back
    ( SpeciesRecord
      ( NAME("HCN"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("124"	,0.98511	,27.	,231	,231	,TAGS(27001, 27003)),
	 REC("134"	,0.01107	,28.	,232	,232	,TAGS(28002)),
	 REC("125"	,0.0036217	,28.	,233	,233	,TAGS(28003)),
	 REC("224"	,0.00014773545	,28.	,-1	,-1	,TAGS(28004))
	 ) ) );

  // CH3Cl
  // Abundancies: 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("CH3Cl"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("215"	,0.74894	,50.	,241	,241	,TAGS(50007)),
	 REC("217"	,0.23949	,52.	,242	,242	,TAGS(52009))
	 ) ) );

  // H2O2
  // Abundancies: 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("H2O2"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("1661"	,0.994952	,34.	,251	,251	,TAGS(34004))
	 ) ) );

  // C2H2
  // Degrees of freedom: guessed, since it seems to be non linear
  // Abundancies: 1 1
  // Note: C2H2 is in official mytran list (26), but does not 
  //       seem to be included for calculation, as given
  //       by table tag_table in file glob_def.c
  species_data.push_back
    ( SpeciesRecord
      ( NAME("C2H2"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("1221"	,0.97760	,26.	,-1	,261	,TAGS( )),
	 REC("1231"	,0.02197	,27.	,-1	,262	,TAGS( ))
	 ) ) );

  // C2H6
  // Degrees of freedom: guessed, since it seems to be non linear
  // Abundancies: 1
  // Note: C2H6 is in official mytran list (27), but does not 
  //       seem to be included for calculation, as given
  //       by table tag_table in file glob_def.c
  species_data.push_back
    ( SpeciesRecord
      ( NAME("C2H6"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("1221"	,0.97699	,30.	,-1	,271	,TAGS( ))
	 ) ) );

  // PH3
  // Abundancies: 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("PH3"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("1111"	,0.99953283	,34.	,281	,281	,TAGS(34003))
	 ) ) );

  // COF2
  // Abundancies: 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("COF2"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("269"	,0.98654	,66.	,291	,291	,TAGS(66001))
	 ) ) );

  // SF6
  // Degrees of freedom: guessed, since it seems to be non linear
  // Abundancies: 1
  // Note: SF6 is in official mytran list (30), but does not 
  //       seem to be included for calculation, as given
  //       by table tag_table in file glob_def.c
  species_data.push_back
    ( SpeciesRecord
      ( NAME("SF6"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("29"	,0.95018	,146.	,-1	,301	,TAGS( ))
	 ) ) );

  // H2S
  // Abundancies: 1 1 1 2
  species_data.push_back
    ( SpeciesRecord
      ( NAME("H2S"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("121"	,0.94988	,34.	,311	,311	,TAGS(34002)),
	 REC("141"	,0.04214	,36.	,-1	,312	,TAGS( )),
	 REC("131"	,0.007498	,35.	,-1	,313	,TAGS( )),
	 REC("122"	,0.00029991625	,35.	,-1	,-1	,TAGS(35001))
	 ) ) );

  // HCOOH
  // Abundancies: 1 3 3 3
  // Note: the isotope names differ from hitran convention, since the jpl catalogue has 
  //       isotopes HCOOH, HC-13-OOH, DCOOH, HCOOD
  //       hitran convention  --  new convention  -- jpl species
  //              126                 1261             HCOOH
  //               -                  1361             HC-13-OOH
  //               -                  2261             DCOOH
  //               -                  1262             HCOOD
  species_data.push_back
    ( SpeciesRecord
      ( NAME("HCOOH"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("1261"	,0.983898	,46.	,321	,321	,TAGS(46005)),
	 REC("1361"	,0.010913149	,47.	,-1	,-1	,TAGS(47002)),
	 REC("2261"	,0.00014755369	,47.	,-1	,-1	,TAGS(47003)),
	 REC("1262"	,0.00014755369	,47.	,-1	,-1	,TAGS(47004))
	 ) ) );

  // HO2
  // Abundancies: 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("HO2"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("166"	,0.995107	,33.	,331	,331	,TAGS(33001))
	 ) ) );

  // O
  // Abundancies: 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("O"),
	DEGFR(0),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("6"	,0.997628	,16.	,341	,341	,TAGS(16001))
	 ) ) );

  // ClONO2
  // Abundancies: 1 1
  // Note: ClONO2 in hitran is identical to ClNO3 in jpl (according to Johannes Orphal)
  species_data.push_back
    ( SpeciesRecord
      ( NAME("ClONO2"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("5646"	,0.74957	,97.	,351	,351	,TAGS(97002)),
	 REC("7646"	,0.23970	,99.	,352	,352	,TAGS(99001))
	 ) ) );

  // NO+
  // Degrees of freedom: guessed, since it seems to be linear
  // Abundancies: 1
  // Note: NO+ is in official mytran list (36), but does not 
  //       seem to be included for calculation, as given
  //       by table tag_table in file glob_def.c
  species_data.push_back
    ( SpeciesRecord
      ( NAME("NO+"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("46"	,0.993974	,30.	,-1	,361	,TAGS(30011))
	 ) ) );

  // OClO
  // Abundancies: 2 2
  species_data.push_back
    ( SpeciesRecord
      ( NAME("OClO"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("656"	,0.75509223	,67.	,431	,-1	,TAGS(67001)),
	 REC("676"	,0.24490632	,69.	,432	,-1	,TAGS(69001))
	 ) ) );

  // BrO
  // Abundancies: 2 2
  species_data.push_back
    ( SpeciesRecord
      ( NAME("BrO"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("96"	,0.50582466	,95.	,461	,-1	,TAGS(95001)),
	 REC("16"	,0.49431069	,97.	,462	,-1	,TAGS(97001))
	 ) ) );

  // H2SO4
  // Abundancies: 2
  species_data.push_back
    ( SpeciesRecord
      ( NAME("H2SO4"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("126"	,0.95060479	,98.	,481	,-1	,TAGS(98001))
	 ) ) );

  // Cl2O2
  // Abundancies: 2 2
  // Note: refered to as Cl2O2 in mytran catalogue, in jpl cat: ClOOCl
  species_data.push_back
    ( SpeciesRecord
      ( NAME("Cl2O2"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("565"	,0.57016427	,102.	,491	,-1	,TAGS(102001)),
	 REC("765"	,0.36982818	,104.	,492	,-1	,TAGS(104001))
	 ) ) );


}
