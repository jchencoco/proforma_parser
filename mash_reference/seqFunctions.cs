using System;
using System.Collections.Generic;
using System.Data;
using System.Drawing;
using System.Text;
using System.Windows.Forms;
using System.Threading;
using System.Collections;
using System.Text.RegularExpressions;

namespace MashMain
{
    public class seqFunctions
    {        
        public struct Mass
        {
            public string sym;
            public string sym_short;
            public string name;
            public string formula;
            public double mono_mass;
            public double avg_mass;

        }

        public static readonly Mass[] AA_RES; // amino acid residue

        public static readonly Mass[] ELEMENTS; // chemical elements

        public static readonly Mass[] COMPOUNDS; // chemical compounds

        static seqFunctions()
        {
            #region set amino acid mass info
            AA_RES = new Mass[20]; // create array

            AA_RES[0].sym = "Ala";
            AA_RES[0].sym_short = "A";
            AA_RES[0].name = "Alanine";
            AA_RES[0].formula = "C3H5N1O1"; //residual formula
            AA_RES[0].mono_mass = 71.03711; //residual mass
            AA_RES[0].avg_mass = 71.0788;

            AA_RES[1].sym = "Arg";
            AA_RES[1].sym_short = "R";
            AA_RES[1].name = "Arginine";
            AA_RES[1].formula = "C6H12N4O1";
            AA_RES[1].mono_mass = 156.10111;
            AA_RES[1].avg_mass = 156.1876;

            AA_RES[2].sym = "Val";
            AA_RES[2].sym_short = "V";
            AA_RES[2].name = "Valine";
            AA_RES[2].formula = "C5H9N1O1";
            AA_RES[2].mono_mass = 99.06841;
            AA_RES[2].avg_mass = 99.1326;

            AA_RES[3].sym = "Asn";
            AA_RES[3].sym_short = "N";
            AA_RES[3].name = "Asparagine";
            AA_RES[3].formula = "C4H6N2O2";
            AA_RES[3].mono_mass = 114.04293;
            AA_RES[3].avg_mass = 114.1039;

            AA_RES[4].sym = "Asp";
            AA_RES[4].sym_short = "D";
            AA_RES[4].name = "Aspartic Acid";
            AA_RES[4].formula = "C4H5N1O3";
            AA_RES[4].mono_mass = 115.02694;
            AA_RES[4].avg_mass = 115.0886;

            AA_RES[5].sym = "Cys";
            AA_RES[5].sym_short = "C";
            AA_RES[5].name = "Cysteine";
            AA_RES[5].formula = "C3H5N1O1S1";
            AA_RES[5].mono_mass = 103.00919;
            AA_RES[5].avg_mass = 103.1448;

            AA_RES[6].sym = "Glu";
            AA_RES[6].sym_short = "E";
            AA_RES[6].name = "Glutamic Acid";
            AA_RES[6].formula = "C5H7N1O3";
            AA_RES[6].mono_mass = 129.04259;
            AA_RES[6].avg_mass = 129.1155;

            AA_RES[7].sym = "Gln";
            AA_RES[7].sym_short = "Q";
            AA_RES[7].name = "Glutamine";
            AA_RES[7].formula = "C5H8N2O2";
            AA_RES[7].mono_mass = 128.05858;
            AA_RES[7].avg_mass = 128.1308;

            AA_RES[8].sym = "Gly";
            AA_RES[8].sym_short = "G";
            AA_RES[8].name = "Glycine";
            AA_RES[8].formula = "C2H3N1O1";
            AA_RES[8].mono_mass = 57.02146;
            AA_RES[8].avg_mass = 57.0520;

            AA_RES[9].sym = "His";
            AA_RES[9].sym_short = "H";
            AA_RES[9].name = "Histidine";
            AA_RES[9].formula = "C6H7N3O1";
            AA_RES[9].mono_mass = 137.05891;
            AA_RES[9].avg_mass = 137.1412;

            AA_RES[10].sym = "Ile";
            AA_RES[10].sym_short = "I";
            AA_RES[10].name = "Isoleucine";
            AA_RES[10].formula = "C6H11N1O1";
            AA_RES[10].mono_mass = 113.08406;
            AA_RES[10].avg_mass = 113.1595;

            AA_RES[11].sym = "Leu";
            AA_RES[11].sym_short = "L";
            AA_RES[11].name = "Leucine";
            AA_RES[11].formula = "C6H11N1O1";
            AA_RES[11].mono_mass = 113.08406;
            AA_RES[11].avg_mass = 113.1595;

            AA_RES[12].sym = "Lys";
            AA_RES[12].sym_short = "K";
            AA_RES[12].name = "Lysine";
            AA_RES[12].formula = "C6H12N2O1";
            AA_RES[12].mono_mass = 128.09496;
            AA_RES[12].avg_mass = 128.1742;

            AA_RES[13].sym = "Met";
            AA_RES[13].sym_short = "M";
            AA_RES[13].name = "Methionine";
            AA_RES[13].formula = "C5H9N1O1S1";
            AA_RES[13].mono_mass = 131.04049;
            AA_RES[13].avg_mass = 131.1986;

            AA_RES[14].sym = "Phe";
            AA_RES[14].sym_short = "F";
            AA_RES[14].name = "Phenylalanine";
            AA_RES[14].formula = "C9H9N1O1";
            AA_RES[14].mono_mass = 147.06841;
            AA_RES[14].avg_mass = 147.1766;

            AA_RES[15].sym = "Pro";
            AA_RES[15].sym_short = "P";
            AA_RES[15].name = "Proline";
            AA_RES[15].formula = "C5H7N1O1";
            AA_RES[15].mono_mass = 97.05276;
            AA_RES[15].avg_mass = 97.1167;

            AA_RES[16].sym = "Ser";
            AA_RES[16].sym_short = "S";
            AA_RES[16].name = "Serine";
            AA_RES[16].formula = "C3H5N1O2";
            AA_RES[16].mono_mass = 87.03203;
            AA_RES[16].avg_mass = 87.0782;

            AA_RES[17].sym = "Thr";
            AA_RES[17].sym_short = "T";
            AA_RES[17].name = "Threonine";
            AA_RES[17].formula = "C4H7N1O2";
            AA_RES[17].mono_mass = 101.04768;
            AA_RES[17].avg_mass = 101.1051;

            AA_RES[18].sym = "Trp";
            AA_RES[18].sym_short = "W";
            AA_RES[18].name = "Tryptophan";
            AA_RES[18].formula = "C11H10N2O1";
            AA_RES[18].mono_mass = 186.07931;
            AA_RES[18].avg_mass = 186.2133;

            AA_RES[19].sym = "Tyr";
            AA_RES[19].sym_short = "Y";
            AA_RES[19].name = "Tyrosine";
            AA_RES[19].formula = "C9H9N1O2";
            AA_RES[19].mono_mass = 163.06333;
            AA_RES[19].avg_mass = 163.1760;
            #endregion

            #region set chemical element mass info
            ELEMENTS = new Mass[11]; // create array

            ELEMENTS[0].sym = "C";
            ELEMENTS[0].sym_short = ELEMENTS[0].sym;
            ELEMENTS[0].name = "Carbon";
            ELEMENTS[0].mono_mass = 12;
            ELEMENTS[0].avg_mass = 12.011;

            ELEMENTS[1].sym = "H";
            ELEMENTS[1].sym_short = ELEMENTS[1].sym;
            ELEMENTS[1].name = "Hydrogen";
            ELEMENTS[1].mono_mass = 1.0078250;
            ELEMENTS[1].avg_mass = 1.00794;

            ELEMENTS[2].sym = "N";
            ELEMENTS[2].sym_short = ELEMENTS[2].sym;
            ELEMENTS[2].name = "Nitrogen";
            ELEMENTS[2].mono_mass = 14.0030740;
            ELEMENTS[2].avg_mass = 14.00674;

            ELEMENTS[3].sym = "O";
            ELEMENTS[3].sym_short = ELEMENTS[3].sym;
            ELEMENTS[3].name = "Oxygen";
            ELEMENTS[3].mono_mass = 15.9949146;
            ELEMENTS[3].avg_mass = 15.9994;

            ELEMENTS[4].sym = "F";
            ELEMENTS[4].sym_short = ELEMENTS[4].sym;
            ELEMENTS[4].name = "Fluorine";
            ELEMENTS[4].mono_mass = 18.9984033;
            ELEMENTS[4].avg_mass = 18.9984;

            ELEMENTS[5].sym = "P";
            ELEMENTS[5].sym_short = ELEMENTS[5].sym;
            ELEMENTS[5].name = "Phosphorous";
            ELEMENTS[5].mono_mass = 30.9737634;
            ELEMENTS[5].avg_mass = 30.97376;

            ELEMENTS[6].sym = "S";
            ELEMENTS[6].sym_short = ELEMENTS[6].sym;
            ELEMENTS[6].name = "Sulfer";
            ELEMENTS[6].mono_mass = 31.9720718;
            ELEMENTS[6].avg_mass = 32.066;

            ELEMENTS[7].sym = "Cl";
            ELEMENTS[7].sym_short = ELEMENTS[7].sym;
            ELEMENTS[7].name = "Chlorine";
            ELEMENTS[7].mono_mass = 34.9688527;
            ELEMENTS[7].avg_mass = 35.4527;

            ELEMENTS[8].sym = "Br";
            ELEMENTS[8].sym_short = ELEMENTS[8].sym;
            ELEMENTS[8].name = "Bromine";
            ELEMENTS[8].mono_mass = 78.9183361;
            ELEMENTS[8].avg_mass = 79.904;

            ELEMENTS[9].sym = "proton";
            ELEMENTS[9].sym_short = ELEMENTS[9].sym;
            ELEMENTS[9].name = "proton";
            ELEMENTS[9].mono_mass = 1.0078250 - 0.0005486; // H - electron
            ELEMENTS[9].avg_mass = 1.00794 - 0.0005486;

            ELEMENTS[10].sym = "el";
            ELEMENTS[10].sym_short = ELEMENTS[10].sym;
            ELEMENTS[10].name = "electron";
            ELEMENTS[10].mono_mass = 0.0005486;
            ELEMENTS[10].avg_mass = 0.0005486;
            #endregion

            #region set chemical compound mass info
            COMPOUNDS = new Mass[4]; // create array

            COMPOUNDS[0].sym = "NH3";
            COMPOUNDS[0].sym_short = COMPOUNDS[0].sym;
            COMPOUNDS[0].name = "Ammonia";
            COMPOUNDS[0].mono_mass = 3 * 1.0078250 + 14.0030740;
            COMPOUNDS[0].avg_mass = 3 * 1.00794 + 14.00674;

            COMPOUNDS[1].sym = "H2O";
            COMPOUNDS[1].sym_short = COMPOUNDS[1].sym;
            COMPOUNDS[1].name = "Water";
            COMPOUNDS[1].mono_mass = 2 * 1.0078250 + 15.9949146;
            COMPOUNDS[1].avg_mass = 2 * 1.00794 + 15.9994;

            COMPOUNDS[2].sym = "HPO3";
            COMPOUNDS[2].sym_short = COMPOUNDS[2].sym;
            COMPOUNDS[2].name = "Phosphorylation";
            COMPOUNDS[2].mono_mass = (1 * 1.0078250) +(1* 30.9737634)+(3*15.9949146);
            COMPOUNDS[2].avg_mass = (1 * 1.00794) +(1*30.97376)+ (3*15.9994);

            COMPOUNDS[3].sym = "H3PO4";
            COMPOUNDS[3].sym_short = COMPOUNDS[2].sym;
            COMPOUNDS[3].name = "Phosphorylation";
            COMPOUNDS[3].mono_mass = (3* 1.0078250) + (1 * 30.9737634) + (4 * 15.9949146);
            COMPOUNDS[3].avg_mass = (3 * 1.00794) + (1 * 30.97376) + (4 * 15.9994);
            #endregion

        }

        public double[] seqTo_N_ion(ref seqTool.ProteinSeq proteinSeq, seqTool.IonType ionType)
        {
            #region conversion factors
            // write algorithm for b ion
            // add uniform shift to change into c or a ion
            double uniform_massShift = 0;

            switch(ionType)
            {
                case seqTool.IonType.a:
                    // need to find conv factors
                    break;
                case seqTool.IonType.b:
                    // do nothing since algorithm is for b
                    break;
                case seqTool.IonType.c:
                    uniform_massShift += getCmpdMass("NH3", proteinSeq.param.massType); // NH3
                    break;
            }

            // add correction for charged ions, is this correct for z>1?
            // effective this is just adding z*proton
            uniform_massShift += (double)proteinSeq.param.charge * getEleMass("proton", proteinSeq.param.massType);
            #endregion

            // convert protein sequence into N ion masses
            double[] N_ion = new double[proteinSeq.sequence.Length];
            N_ion[0] = getAaMass(proteinSeq.sequence[0],proteinSeq.param.massType) + uniform_massShift + proteinSeq.ptm.massChange[0];
            for (int i = 1; i < N_ion.Length; i++)
            {
                N_ion[i] = N_ion[i - 1] + getAaMass(proteinSeq.sequence[i], proteinSeq.param.massType) + proteinSeq.ptm.massChange[i];
            }

            return N_ion;
        }

        public double[] seqTo_C_ion(ref seqTool.ProteinSeq proteinSeq, seqTool.IonType ionType)
        {
            #region conversion factors
            // write algorithm for y ion
            // add uniform shift to change into z or x ion
            double uniform_massShift = 0;

            switch (ionType)
            {
                case seqTool.IonType.x:
                    // need to find conv factors
                    break;
                case seqTool.IonType.y:
                    // do nothing since alg. is written for y ions
                    break;
                case seqTool.IonType.z:
                    uniform_massShift -= 2 * getEleMass("H", proteinSeq.param.massType) + getEleMass("N", proteinSeq.param.massType); // -NH2
                    break;
            }

            // add correction for charged ions, is this correct for z>1?
            // effective this is just adding z*proton
            uniform_massShift += (double)proteinSeq.param.charge * getEleMass("proton", proteinSeq.param.massType);
            #endregion

            // convert protein sequence into y ion masses
            int max = proteinSeq.sequence.Length;

            double[] C_ion = new double[max];
            double water = getCmpdMass("H2O", proteinSeq.param.massType);

            C_ion[0] = getAaMass(proteinSeq.sequence[max - 1], proteinSeq.param.massType) + water + uniform_massShift + proteinSeq.ptm.massChange[max - 1];

            for (int i = 1; i < C_ion.Length; i++)
            {
                int curIndex = max - i - 1;
                C_ion[i] = C_ion[i - 1] + getAaMass(proteinSeq.sequence[curIndex], proteinSeq.param.massType) + proteinSeq.ptm.massChange[curIndex];
            }
            return C_ion;
        }

        public List<seqTool.ProteinSeq.Peptide> computePeptideMass(seqTool.ProteinSeq proteinSeq)
        {
            // compute peptide mass with ptm
            // note: peptide mass = sum(aa_res mass) + ptm(aa_res) + H2O <- H at N-terminal, OH at C-terminal
              seqTool.ProteinSeq.Peptide[] pArray = proteinSeq.peptide.ToArray();

            for (int i=0; i<pArray.Length; i++)
            {
                int absInd;
                for (int j=0; j<pArray[i].sequence.Length; j++)
                {
                    absInd = j + pArray[i].firstInd - 1;
                    pArray[i].mass += getAaMass(pArray[i].sequence[j],proteinSeq.param.massType) + proteinSeq.ptm.massChange[absInd];
                }
                pArray[i].mass += getCmpdMass("H2O", proteinSeq.param.massType);
                pArray[i].mass += (double)proteinSeq.param.charge * getEleMass("proton", proteinSeq.param.massType); // for [M+H]+
            }

            List<seqTool.ProteinSeq.Peptide> pList = new List<seqTool.ProteinSeq.Peptide>();
            pList.AddRange(pArray);
            return pList;
        }

        public string getFormulaSequence(string proteinSeq)
        {
         

            //for (int i = 0; i < proteinSeq.Length; i++)
            //{
            //    int absInd;

            //}

            return "CH2";
        }

        public void compareMass(ref seqTool.ProteinSeq proteinSeq, seqTool.IonType ionType)
        {
            // compares one theoretial ion type against experimental value



            #region select the correct theoretical ion mass to use for comparison
            double[] massTheo;
            switch(ionType)
            {
                case seqTool.IonType.a:
                    massTheo = proteinSeq.mass.aIon;
                    break;
                case seqTool.IonType.b:
                    massTheo = proteinSeq.mass.bIon;
                    break;
                case seqTool.IonType.c:
                    massTheo = proteinSeq.mass.cIon;
                    break;
                case seqTool.IonType.x:
                    massTheo = proteinSeq.mass.xIon;
                    break;
                case seqTool.IonType.y:
                    massTheo = proteinSeq.mass.yIon;
                    break;
                case seqTool.IonType.z:
                    massTheo = proteinSeq.mass.zIon;
                    break;
                default:
                    massTheo = new double[proteinSeq.sequence.Length];
                    break;
                    // not necessary but C# complains about unassigned variable..
                    // alternatively can make one of the cases the default one, such as zIon
            }
            #endregion



            double[] massExptl = proteinSeq.mass.exprtl;

            string ionOutputPrefix = ionType.ToString();
            double tolerance;

            // for each experimental mass, loop through ionMass to find the best match
            for (int i = 0; i < massExptl.Length; i++)
            {
                #region see if experimental mass has any old matches
                double lowestDiff;

                if (proteinSeq.seqResult[i].massTheoretical == 0) // no previous matches, new comparison
                {
                    proteinSeq.seqResult[i].ionAssignment = "NA";
                    lowestDiff = 10000;
                }
                else // see if this iontype is any better
                {
                    lowestDiff = proteinSeq.seqResult[i].error;
                }
                #endregion

                #region iterate through theoretical ion mass list to find best match for experimental mass
                for (int j = 0; j < massTheo.Length; j++)
                {
                    #region Convert tolerance to Da
                    // Note: comparison is using 2*tolerance due to Abs()


                    // convert ppm to Da
                    if (proteinSeq.param.tolType == seqTool.ToleranceType.ppm)
                    {
                        tolerance = proteinSeq.param.tolerance * massTheo[j] / 1E+6;
                    }
                    else // Da, no conversion needed
                    {
                        tolerance = proteinSeq.param.tolerance;
                    }
                    #endregion

                    double diff = massExptl[i] - massTheo[j];
                    if (Math.Abs(diff) <= tolerance && Math.Abs(diff) < lowestDiff)
                    {
                        proteinSeq.seqResult[i].ionAssignment = ionOutputPrefix + (j + 1).ToString();
                        proteinSeq.seqResult[i].massTheoretical = massTheo[j];
                        proteinSeq.seqResult[i].error = Math.Round(diff, 5);
                        lowestDiff = diff;
                    }

                }
                #endregion
            }

        }

        public void compareMass(ref seqTool.ProteinSeq proteinSeq, double[] massTheo, string outputPrefix)
        {
            // compares one theoretial ion type against experimental value

            double[] massExptl = proteinSeq.mass.exprtl;
            double tolerance;

            // for each experimental mass, loop through ionMass to find the best match
            for (int i = 0; i < massExptl.Length; i++)
            {
                #region see if experimental mass has any old matches
                double lowestDiff;

                if (proteinSeq.seqResult[i].massTheoretical == 0) // no previous matches, new comparison
                {
                    proteinSeq.seqResult[i].ionAssignment = "NA";
                    lowestDiff = 10000;
                }
                else // see if this iontype is any better
                {
                    lowestDiff = proteinSeq.seqResult[i].error;
                }
                #endregion

                #region iterate through theoretical ion mass list to find best match for experimental mass
                for (int j = 0; j < massTheo.Length; j++)
                {
                    #region Convert tolerance to Da
                    // Note: comparison is using 2*tolerance due to Abs()


                    // convert ppm to Da
                    if (proteinSeq.param.tolType == seqTool.ToleranceType.ppm)
                    {
                        tolerance = proteinSeq.param.tolerance * massTheo[j] / 1E+6;
                    }
                    else // Da, no conversion needed
                    {
                        tolerance = proteinSeq.param.tolerance;
                    }
                    #endregion

                    double diff = massExptl[i] - massTheo[j];
                    if (Math.Abs(diff) <= tolerance && Math.Abs(diff) < lowestDiff)
                    {
                        proteinSeq.seqResult[i].ionAssignment = outputPrefix + (j + 1).ToString();
                        proteinSeq.seqResult[i].massTheoretical = massTheo[j];
                        proteinSeq.seqResult[i].error = Math.Round(diff, 5);
                        lowestDiff = diff;
                    }

                }
                #endregion
            }

        }

        public double getAaMass(char aa, seqTool.MassType massType)
        {
            for (int i = 0; i < AA_RES.Length; i++)
                if (string.Equals(AA_RES[i].sym_short, aa.ToString()))
                    if (massType == seqTool.MassType.monoisotopic)
                        return AA_RES[i].mono_mass;
                    else
                        return AA_RES[i].avg_mass;

            // return negative value if not found?
            return -1;
        }

        public double getEleMass(string sym, seqTool.MassType massType)
        {
            for (int i = 0; i < ELEMENTS.Length; i++)
                if (string.Equals(ELEMENTS[i].sym, sym))
                    if (massType == seqTool.MassType.monoisotopic)
                        return ELEMENTS[i].mono_mass;
                    else
                        return ELEMENTS[i].avg_mass;

            // return negative value if not found?
            return -1;
        }

        public double getCmpdMass(string sym, seqTool.MassType massType)
        {
            for (int i = 0; i < COMPOUNDS.Length; i++)
                if (string.Equals(COMPOUNDS[i].sym, sym))
                    if (massType == seqTool.MassType.monoisotopic)
                        return COMPOUNDS[i].mono_mass;
                    else
                        return COMPOUNDS[i].avg_mass;

            // return negative value if not found?
            return -1;
        }

        public double[] massShift(double shiftAmt, double[] massArray)
        {
            // usually shiftAmt is negative. eg. water and ammonia loss

            for (int i = 0; i < massArray.Length; i++)
            {
                massArray[i] += shiftAmt;
            }

            return massArray;
        }

        public List<seqTool.ProteinSeq.Peptide> cleaveWithEnzyme(string peptideSeq, string enzyme, int numMissed)
        {
            // find cleavage index
            int[] validInd = findCleaveInd(peptideSeq, enzyme);

            // cleave on C-terminal side
            List<seqTool.ProteinSeq.Peptide> cleavedPeptide = cleavePeptide(peptideSeq, validInd);

            // generate missed cleavage fragments from fully cleaved fragments
            cleavedPeptide = mcPeptide(cleavedPeptide, numMissed);

            return cleavedPeptide;
        }

        private int[] findCleaveInd(string peptideSeq, string enzyme)
        {
            // index for first aa_res is 0
            Regex cleaveSite;
            MatchCollection matches;
            ArrayList validInd;
            bool isValid;
            int[] sortedInd;

            switch (enzyme)
            {
                case "trypsin":
                    cleaveSite = new Regex("[KR]");
                    matches = cleaveSite.Matches(peptideSeq);
                    validInd = new ArrayList();

                    foreach (Match match in matches)
                    {
                        isValid = false;

                        int nextCind = match.Index + 1;

                        // Exception: if P is C-term to K or R
                        if (nextCind >= peptideSeq.Length)
                            isValid = true;
                        else if (peptideSeq[nextCind] != 'P')
                            isValid = true;

                        if (isValid)
                            validInd.Add(match.Index);
                    }

                    validInd.Sort();
                    sortedInd = (int[]) validInd.ToArray(typeof(int));
                    return sortedInd;
                case "lys c":
                    cleaveSite = new Regex("[K]");
                    matches = cleaveSite.Matches(peptideSeq);
                    validInd = new ArrayList();

                    // Exception: none
                    foreach (Match match in matches)
                            validInd.Add(match.Index);
                    
                    validInd.Sort();
                    sortedInd = (int[])validInd.ToArray(typeof(int));
                    return sortedInd;

                case "asp n":
                    cleaveSite = new Regex("[D]");
                    matches = cleaveSite.Matches(peptideSeq);
                    validInd = new ArrayList();

                    // Exception: none
                    foreach (Match match in matches)
                        validInd.Add(match.Index - 1);

                    validInd.Sort();
                    sortedInd = (int[])validInd.ToArray(typeof(int));
                    return sortedInd;

                case "lys n":
                    cleaveSite = new Regex("[K]");
                    matches = cleaveSite.Matches(peptideSeq);
                    validInd = new ArrayList();

                    // Exception: none
                    foreach (Match match in matches)
                        validInd.Add(match.Index - 1);

                    validInd.Sort();
                    sortedInd = (int[])validInd.ToArray(typeof(int));
                    return sortedInd;

                case "cnbr":
                    cleaveSite = new Regex("[M]");
                    matches = cleaveSite.Matches(peptideSeq);
                    validInd = new ArrayList();

                    // Exception: none
                    foreach (Match match in matches)
                        validInd.Add(match.Index);
                 
                    validInd.Sort();
                    sortedInd = (int[])validInd.ToArray(typeof(int));
                    return sortedInd;
                case "glu c (bicarbonate)":
                    cleaveSite = new Regex("[E]");
                    matches = cleaveSite.Matches(peptideSeq);
                    validInd = new ArrayList();

                    foreach (Match match in matches)
                    {
                        isValid = true;

                        int nextCind = match.Index + 1;
                        // Exception: if P is C-term to E, or if E is C-term to E

                        if (nextCind >= peptideSeq.Length)
                            isValid = true;
                        else if (peptideSeq[nextCind] == 'P')
                            isValid = false;
                        else if (peptideSeq[nextCind] == 'E')
                            isValid = false;

                        if (isValid)
                            validInd.Add(match.Index);
                    }

                    validInd.Sort();
                    sortedInd = (int[])validInd.ToArray(typeof(int));
                    return sortedInd;

                case "arg c":
                    cleaveSite = new Regex("[R]");
                    matches = cleaveSite.Matches(peptideSeq);
                    validInd = new ArrayList();

                    foreach (Match match in matches)
                    {
                        isValid = true;

                        int nextCind = match.Index +1;
                        // Exception if if P is C-term to R
                        if (nextCind >= peptideSeq.Length)
                                isValid = true;
                        else if (peptideSeq[nextCind] == 'P')
                                isValid = false;
                        if (isValid)
                                validInd.Add(match.Index);

                        }
                    validInd.Sort();
                    sortedInd = (int[])validInd.ToArray(typeof(int));
                    return sortedInd;

                case "pepsin (ph 1.3)":
                    cleaveSite = new Regex("[FL]");
                    matches = cleaveSite.Matches(peptideSeq);
                    validInd = new ArrayList();

                    foreach (Match match in matches)
                        validInd.Add(match.Index);
                    //Exception None

                    validInd.Sort();
                    sortedInd = (int[])validInd.ToArray(typeof(int));
                    return sortedInd;

                case "pepsin (ph > 2)":
                    cleaveSite = new Regex("[FLWYAEQ]");
                    matches = cleaveSite.Matches(peptideSeq);
                    validInd = new ArrayList();

                    foreach (Match match in matches)
                        validInd.Add(match.Index);
                    //Exception None

                    validInd.Sort();
                    sortedInd = (int[])validInd.ToArray(typeof(int));
                    return sortedInd;

                case "proteinase k":
                    cleaveSite = new Regex("[AFYWLIV]");
                    matches = cleaveSite.Matches(peptideSeq);
                    validInd = new ArrayList();

                    foreach (Match match in matches)
                        validInd.Add(match.Index);
                    //Exception None

                    validInd.Sort();
                    sortedInd = (int[])validInd.ToArray(typeof(int));
                    return sortedInd;

                case "thermolysin":
                    cleaveSite = new Regex("[AFILMV]");
                    matches = cleaveSite.Matches(peptideSeq);
                    validInd = new ArrayList();

                    foreach (Match match in matches)
                    {
                        isValid = true;
                        
                        int nextNind = match.Index + 1;
                        // Exception: if D or E is N-term to A, F, I, L, M, V
                        if (nextNind == peptideSeq.Length)
                           isValid = true;
                        if (peptideSeq[nextNind - 2] == 'D')
                            isValid = false;
                        else if (peptideSeq[nextNind - 2] == 'E')
                            isValid = false;

                        if (isValid)
                            validInd.Add(match.Index - 1);
                      
                    }

                    validInd.Sort();
                    sortedInd = (int[])validInd.ToArray(typeof(int));
                    return sortedInd;                    

                case "glu c (phosphate)":
                    cleaveSite = new Regex("[DE]");
                    matches = cleaveSite.Matches(peptideSeq);
                    validInd = new ArrayList();

                    foreach (Match match in matches)
                    {
                        isValid = true;

                        int nextCind = match.Index + 1;
                        // Exception: if P is C-term to D or E, or if E is C-term to D or E
                        if (nextCind >= peptideSeq.Length)
                            isValid = true;
                        else if (peptideSeq[nextCind] == 'P')
                            isValid = false;
                        else if (peptideSeq[nextCind] == 'E')
                            isValid = false;

                        if (isValid)
                            validInd.Add(match.Index);
                    }

                    validInd.Sort();
                    sortedInd = (int[])validInd.ToArray(typeof(int));
                    return sortedInd;
                default:
                    break;
            }

            return null;
        }

        private List<seqTool.ProteinSeq.Peptide> cleavePeptide(string peptideSeq, int[] ind)
        {
            int leftInd = 0;
            int len;

            List<seqTool.ProteinSeq.Peptide> cleavedPeptide = new List<seqTool.ProteinSeq.Peptide>();

            for (int i = 0; i < ind.Length; i++)
            {
                len = ind[i] - leftInd + 1;
                seqTool.ProteinSeq.Peptide temp = new seqTool.ProteinSeq.Peptide();
                temp.sequence = peptideSeq.Substring(leftInd, len);
                temp.firstInd = leftInd + 1;
                temp.lastInd = leftInd + len;
                temp.mc = 0;
                cleavedPeptide.Add(temp);
                leftInd = leftInd + len;
            }

            // last cleaved piece
            if (leftInd < peptideSeq.Length-1)
            {
                len = peptideSeq.Length - leftInd;
                seqTool.ProteinSeq.Peptide temp = new seqTool.ProteinSeq.Peptide();
                temp.sequence = peptideSeq.Substring(leftInd, len);
                temp.firstInd = leftInd + 1;
                temp.lastInd = leftInd + len;
                temp.mc = 0;
                cleavedPeptide.Add(temp);
            }

            return cleavedPeptide;
        }

        private List<seqTool.ProteinSeq.Peptide> mcPeptide(List<seqTool.ProteinSeq.Peptide> cleavedPeptide, int numMissed)
        {
            int maxNumCleavages = cleavedPeptide.Count - 1;

            // cannot miss more than max # of cleavages possible
            if (numMissed > maxNumCleavages)
                numMissed = maxNumCleavages;

            seqTool.ProteinSeq.Peptide[] prevLvl = cleavedPeptide.ToArray();

            for (int j = 0; j < numMissed; j++) // each mc level
            {
                List<seqTool.ProteinSeq.Peptide> nextLvl = new List<seqTool.ProteinSeq.Peptide>();

                for (int k = 0; k < prevLvl.Length - 1; k++)
                {
                    seqTool.ProteinSeq.Peptide CtermFrag = cleavedPeptide[j + k+1];
                    nextLvl.Add(mergePeptide(prevLvl[k],CtermFrag));
                }

                cleavedPeptide.AddRange(nextLvl);
                prevLvl = nextLvl.ToArray();
            }

            return cleavedPeptide;
        }

        private seqTool.ProteinSeq.Peptide mergePeptide(seqTool.ProteinSeq.Peptide a, seqTool.ProteinSeq.Peptide b)
        {
            seqTool.ProteinSeq.Peptide c = new seqTool.ProteinSeq.Peptide();

            c.sequence = a.sequence + b.sequence;
            c.firstInd = a.firstInd;
            c.lastInd = b.lastInd;
            c.mc = a.mc + 1;

            return c;
        }

        public seqTool.ProteinSeq.PeptideOutput[] comparePeptideMass(seqTool.ProteinSeq proteinSeq)
        {

            double[] massExptl = proteinSeq.mass.exprtl;
            double tolerance;

            seqTool.ProteinSeq.PeptideOutput[] result = new seqTool.ProteinSeq.PeptideOutput[massExptl.Length];

            #region compare masses
            for (int i = 0; i < massExptl.Length; i++)
            {
                #region initialize result
                result[i].sequence = "NA";
                result[i].position = "NA";
                result[i].massExperimental = massExptl[i];
                result[i].massTheoretical = 0;
                result[i].errorDa = 0;
                result[i].errorPPM = 0;
                result[i].mc = 0;
                #endregion

                double lowestDiff = 100000; // arbitrarily high #

                foreach (seqTool.ProteinSeq.Peptide p in proteinSeq.peptide)
                {
                    #region Convert tolerance to Da
                    // Note: comparison is using 2*tolerance due to Abs()
                    // convert ppm to Da
                    if (proteinSeq.param.tolType == seqTool.ToleranceType.ppm)
                        tolerance = proteinSeq.param.tolerance * p.mass / 1E+6;
                    else // Da, no conversion needed
                        tolerance = proteinSeq.param.tolerance;
                    #endregion

                    double diff = massExptl[i] - p.mass;
                    if (Math.Abs(diff) <= tolerance && Math.Abs(diff) < lowestDiff)
                    {
                        result[i].mc = p.mc;
                        result[i].position = p.firstInd.ToString() + "-" + p.lastInd.ToString();
                        result[i].sequence = p.sequence;
                        result[i].massTheoretical = p.mass;
                        result[i].errorDa = Math.Round(diff,5);
                        result[i].errorPPM = Math.Round(diff / p.mass * 1E+6, 5);
                        lowestDiff = diff;
                    }
                }
            }
            #endregion

            return result;
        }
    }
}
