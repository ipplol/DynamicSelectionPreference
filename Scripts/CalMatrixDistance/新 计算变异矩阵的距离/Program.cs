using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Threading;

namespace 新_计算变异矩阵的距离
{
    public class VOC
    {
        public string VOCname;
        public List<string> Weeklist = new List<string>();
        public List<string> WeekCases = new List<string>();
        public List<string> WeekAccumulateCases = new List<string>();
        public List<string> WeekInfectionProp = new List<string>();

        public List<string> SpikeMutList = new List<string>();
        public Dictionary<string, double> DateMutFreqDic = new Dictionary<string, double>();//key = Date + "_" + Mut

        public Dictionary<string, double> Spike_CellEntry = new Dictionary<string, double>();
        public Dictionary<string, double> RBD_ACE2Binding = new Dictionary<string, double>();
        public Dictionary<string, double> RBD_EscapeScore = new Dictionary<string, double>();
    }
    internal class Program
    {
        static string[] VOClist = "WT,Alpha,Delta,BA.1,BA.2,BA.4/5,BQ.1,XBB".Split(',');
        static Dictionary<string, VOC> VOCDic = new Dictionary<string, VOC>();
        static Dictionary<string, double> AntigenicDistance = new Dictionary<string, double>();

        static List<string> PreOmicronEpitope = new List<string>();
        static List<string> OmicronEpitope = new List<string>();
        static Dictionary<int, List<double>> PreOmicronEpitopeDic = new Dictionary<int, List<double>>();//位置与表位信息
        static Dictionary<int, List<double>> OmicronEpitopeDic = new Dictionary<int, List<double>>();
        static Dictionary<string, List<double>> VOCPressureDic = new Dictionary<string, List<double>>();

        static void ReadinDMS()
        {
            int i, j, k;
            string line;
            for(i=0;i<VOClist.Count();i++)
            {
                VOC newv = new VOC();
                newv.VOCname = VOClist[i];
                VOCDic.Add(VOClist[i], newv);
            }

            StreamReader read = new StreamReader("//NAS8500/g/VariationMutation/AiDMS/Spike0417/Data/DMS/Spike_CellEntry.txt");
            line = read.ReadLine(); line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                if (VOCDic.ContainsKey(line1[0]))
                {
                    if (line1[4]!="")
                        VOCDic[line1[0]].Spike_CellEntry.Add(line1[1] + line1[3], Convert.ToDouble(line1[4]));
                    else
                        VOCDic[line1[0]].Spike_CellEntry.Add(line1[1] + line1[3], 0.0);
                }
                line = read.ReadLine();
            }
            read.Close();

            read = new StreamReader("//NAS8500/g/VariationMutation/AiDMS/Spike0417/Data/DMS/RBD_ACE2.txt");
            line = read.ReadLine(); line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split('\t');
                if (VOCDic.ContainsKey(line1[0]))
                {
                    if(line1[2]!="NA")
                        VOCDic[line1[0]].RBD_ACE2Binding.Add(line1[1].Substring(1,4), Convert.ToDouble(line1[2]));
                    else
                        VOCDic[line1[0]].RBD_ACE2Binding.Add(line1[1].Substring(1, 4), 0.0);
                }
                line = read.ReadLine();
            }
            read.Close();

            read = new StreamReader("//NAS8500/g/VariationMutation/AiDMS/Spike0417/Data/DMS/EscapeScore_Merged.txt");
            line = read.ReadLine(); line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split('\t');
                if (VOCDic.ContainsKey(line1[0]))
                {
                    VOCDic[line1[0]].RBD_EscapeScore.Add(line1[1].Substring(1,4), Convert.ToDouble(line1[2]));
                }
                line = read.ReadLine();
            }
            read.Close();

            read = new StreamReader("//NAS8500/g/VariationMutation/AiDMS/Spike0417/Data/DMS/EVEScape_DisAcc_sigmoid.txt");
            line = read.ReadLine(); line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split('\t');
                AntigenicDistance.Add(line1[0] + line1[2], Convert.ToDouble(line1[3]));
                line = read.ReadLine();
            }
            read.Close();
            return;
        }
        static void ReadinMatrix()
        {
            string[] filePaths = Directory.GetFiles("//NAS8500/g/VariationMutation/AiDMS/Spike0417/WeeklyMatrix/Matrix", "*.tsv", SearchOption.TopDirectoryOnly);
            int i, j, k;
            for(k=0;k<filePaths.Count();k++)
            {
                StreamReader read = new StreamReader(filePaths[k]);
                string line = read.ReadLine();
                string[] title = line.Split('\t');
                if (VOCDic.ContainsKey(title[0]))
                    for (i = 6; i < title.Count(); i++)
                        VOCDic[title[0]].SpikeMutList.Add(title[i]);

                    line = read.ReadLine();
                while(line!=null)
                {
                    string[] line1 = line.Split('\t');

                    if (line1[2]=="Y" && VOCDic.ContainsKey(line1[0]))
                    {
                        VOCDic[line1[0]].Weeklist.Add(line1[1]);
                        VOCDic[line1[0]].WeekCases.Add(line1[3]);
                        VOCDic[line1[0]].WeekAccumulateCases.Add(line1[4]);
                        VOCDic[line1[0]].WeekInfectionProp.Add(line1[5]);

                        for (i = 6; i < title.Count(); i++)
                            VOCDic[title[0]].DateMutFreqDic.Add(line1[1] + "_" + title[i], Convert.ToDouble(line1[i]));
                    }
                    line = read.ReadLine();
                }
                read.Close();
            }
            return;
        }
        static void CalculateDistance(string fold)
        {
            int i, j, k;
            StreamWriter writeMerged = new StreamWriter(fold + "/WeeklyDistance.merged.txt");
            writeMerged.WriteLine("VOC\tDate\tWeeklyCases\tInfectionProp\tAntigenicDistanceT1\tSpikeCellEntryT1\tRBDACE2bindingT1\tEscapeScoreT1\tAntigenicDistanceTn-1\tSpikeCellEntryTn-1\tRBDACE2bindingTn-1\tEscapeScoreTn-1");
            for (i=0;i<VOClist.Count();i++)
            {
                StreamWriter write = new StreamWriter(fold + "/WeeklyDistance." + VOClist[i].Replace("/","_") + ".tsv");
                write.WriteLine("VOC\tDate\tWeeklyCases\tInfectionProp\tAntigenicDistanceT1\tSpikeCellEntryT1\tRBDACE2bindingT1\tEscapeScoreT1\tAntigenicDistanceTn-1\tSpikeCellEntryTn-1\tRBDACE2bindingTn-1\tEscapeScoreTn-1");
                for (j = 0; j < VOCDic[VOClist[i]].Weeklist.Count();j++)
                {
                    string output = VOClist[i] + "\t" + VOCDic[VOClist[i]].Weeklist[j] + "\t" + VOCDic[VOClist[i]].WeekCases[j] + "\t" + VOCDic[VOClist[i]].WeekInfectionProp[j];
                    double weeklyAntigenicDistanceT1 = 0, weeklySpikeCellEntryT1 = 0, weeklyRBDACE2bindingT1 = 0, weeklyEscapeScoreT1 = 0;
                    double weeklyAntigenicDistanceTn_1 = 0, weeklySpikeCellEntryTn_1 = 0, weeklyRBDACE2bindingTn_1 = 0, weeklyEscapeScoreTn_1 = 0;
                    for (k=0;k< VOCDic[VOClist[i]].SpikeMutList.Count();k++)
                    {
                        if (VOCDic[VOClist[i]].SpikeMutList[k] == "")
                            continue;

                        string mutk = VOCDic[VOClist[i]].SpikeMutList[k].Substring(1, VOCDic[VOClist[i]].SpikeMutList[k].Length-1);
                        string keyT1 = VOCDic[VOClist[i]].Weeklist[0] + "_" + VOCDic[VOClist[i]].SpikeMutList[k];//以第一个window为参考系
                        double freqT1 = VOCDic[VOClist[i]].DateMutFreqDic[keyT1];
                        string keyTn_1; 
                        if(j>0)
                            keyTn_1 = VOCDic[VOClist[i]].Weeklist[j-1] + "_" + VOCDic[VOClist[i]].SpikeMutList[k];
                        else
                            keyTn_1 = VOCDic[VOClist[i]].Weeklist[0] + "_" + VOCDic[VOClist[i]].SpikeMutList[k];
                        double freqTn_1 = VOCDic[VOClist[i]].DateMutFreqDic[keyTn_1];
                        string keyTn = VOCDic[VOClist[i]].Weeklist[j] + "_" + VOCDic[VOClist[i]].SpikeMutList[k];
                        double freqTn = VOCDic[VOClist[i]].DateMutFreqDic[keyTn];
                        
                        double distanceT1Freq = freqTn - freqT1;
                        double distanceTn_1Freq = freqTn - freqTn_1;

                        //T1
                        if(AntigenicDistance.ContainsKey(mutk))
                            weeklyAntigenicDistanceT1 += AntigenicDistance[mutk] * Math.Abs(distanceT1Freq);
                        if (VOCDic[VOClist[i]].Spike_CellEntry.ContainsKey(mutk))
                            weeklySpikeCellEntryT1 += VOCDic[VOClist[i]].Spike_CellEntry[mutk] * distanceT1Freq;
                        if (VOCDic[VOClist[i]].RBD_ACE2Binding.ContainsKey(mutk))
                            weeklyRBDACE2bindingT1 += VOCDic[VOClist[i]].RBD_ACE2Binding[mutk] * distanceT1Freq;
                        if (VOCDic[VOClist[i]].RBD_EscapeScore.ContainsKey(mutk))
                            weeklyEscapeScoreT1 += VOCDic[VOClist[i]].RBD_EscapeScore[mutk] * distanceT1Freq;
                        //Tn-1
                        if (AntigenicDistance.ContainsKey(mutk))
                            weeklyAntigenicDistanceTn_1 += AntigenicDistance[mutk] * Math.Abs(distanceTn_1Freq);
                        if (VOCDic[VOClist[i]].Spike_CellEntry.ContainsKey(mutk))
                            weeklySpikeCellEntryTn_1 += VOCDic[VOClist[i]].Spike_CellEntry[mutk] * distanceTn_1Freq;
                        if (VOCDic[VOClist[i]].RBD_ACE2Binding.ContainsKey(mutk))
                            weeklyRBDACE2bindingTn_1 += VOCDic[VOClist[i]].RBD_ACE2Binding[mutk] * distanceTn_1Freq;
                        if (VOCDic[VOClist[i]].RBD_EscapeScore.ContainsKey(mutk))
                            weeklyEscapeScoreTn_1 += VOCDic[VOClist[i]].RBD_EscapeScore[mutk] * distanceTn_1Freq;
                    }
                    output += "\t" + Convert.ToString(weeklyAntigenicDistanceT1);
                    output += "\t" + Convert.ToString(weeklySpikeCellEntryT1);
                    output += "\t" + Convert.ToString(weeklyRBDACE2bindingT1);
                    output += "\t" + Convert.ToString(weeklyEscapeScoreT1);
                    output += "\t" + Convert.ToString(weeklyAntigenicDistanceTn_1);
                    output += "\t" + Convert.ToString(weeklySpikeCellEntryTn_1);
                    output += "\t" + Convert.ToString(weeklyRBDACE2bindingTn_1);
                    output += "\t" + Convert.ToString(weeklyEscapeScoreTn_1);
                    write.WriteLine(output);
                    writeMerged.WriteLine(output);
                }
                write.Close();
            }
            writeMerged.Close();
            return;
        }
        static void CalculateAntigenDistanceMatrix(string fold)//计算所有voc所有日期的距离矩阵
        {
            int i, j, k;
            StreamWriter writeMerged = new StreamWriter(fold + "/AllVOCAllWeeklyDistanceMatrix.txt");
            int x, y;
            string output = "AntigenicDistance";
            for (i = 0; i < VOClist.Count(); i++)
                for (j = 0; j < VOCDic[VOClist[i]].Weeklist.Count(); j++)
                    output += "\t" + VOClist[i] + "_" + VOCDic[VOClist[i]].Weeklist[j];
            writeMerged.WriteLine(output);
            
            for (i=0;i<VOClist.Count();i++)
            {
                for (j = 0; j < VOCDic[VOClist[i]].Weeklist.Count();j++)
                {
                    //
                    output = VOClist[i] + "_" + VOCDic[VOClist[i]].Weeklist[j];

                    for (x = 0; x < VOClist.Count(); x++)
                    {
                        for (y = 0; y < VOCDic[VOClist[x]].Weeklist.Count(); y++)
                        {
                            double distance = 0;
                            for (k = 0; k < VOCDic[VOClist[x]].SpikeMutList.Count(); k++)
                            {
                                if (VOCDic[VOClist[x]].SpikeMutList[k] != "")
                                {
                                    double distanceFreq1 = VOCDic[VOClist[i]].DateMutFreqDic[VOCDic[VOClist[i]].Weeklist[j] + "_" + VOCDic[VOClist[i]].SpikeMutList[k]];
                                    double distanceFreq2 = VOCDic[VOClist[x]].DateMutFreqDic[VOCDic[VOClist[x]].Weeklist[y] + "_" + VOCDic[VOClist[x]].SpikeMutList[k]];
                                    string mutk = VOCDic[VOClist[x]].SpikeMutList[k].Substring(1, VOCDic[VOClist[x]].SpikeMutList[k].Length - 1);
                                    if (AntigenicDistance.ContainsKey(mutk))
                                        distance += AntigenicDistance[mutk] * Math.Abs(distanceFreq1 - distanceFreq2);
                                }
                            }
                            output += "\t" + Convert.ToString(distance);
                        }
                    }

                    writeMerged.WriteLine(output);
                }
            }
            writeMerged.Close();
            return;
        }
        static void ReadinEpitope()
        {
            int i, j, k;
            StreamReader read = new StreamReader("//NAS8500/g/VariationMutation/AiDMS/Spike1205/DMSData/EscapeScore_WithNeu_PreOmicron.12.Epitope");
            string line = read.ReadLine();
            string[] title = line.Split('\t');
            for (i = 1; i < title.Count(); i++)
                PreOmicronEpitope.Add(title[i]);
            line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split('\t');
                List<double> newd = new List<double>();
                for (i = 1; i < line1.Count(); i++)
                    newd.Add(Convert.ToDouble(line1[i]));
                PreOmicronEpitopeDic.Add(Convert.ToInt32(line1[0]), newd);
                line = read.ReadLine();
            }
            read.Close();

            read = new StreamReader("//NAS8500/g/VariationMutation/AiDMS/Spike1205/DMSData/EscapeScore_WithNeu_Omicron.12.Epitope");
            line = read.ReadLine();
            title = line.Split('\t');
            for (i = 1; i < title.Count(); i++)
                OmicronEpitope.Add(title[i]);
            line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split('\t');
                List<double> newd = new List<double>();
                for (i = 1; i < line1.Count(); i++)
                    newd.Add(Convert.ToDouble(line1[i]));
                OmicronEpitopeDic.Add(Convert.ToInt32(line1[0]), newd);
                line = read.ReadLine();
            }
            read.Close();

            read = new StreamReader("//NAS8500/g/VariationMutation/AiDMS/Spike0417/WeeklyMatrix/Epitope/AntibodyPressure12Epitope.txt");
            line = read.ReadLine();
            line = read.ReadLine();
            while (line!=null)
            {
                string[] line1 = line.Split('\t');
                List<double> newd = new List<double>();
                for (i = 1; i < line1.Count(); i++)
                    newd.Add(Convert.ToDouble(line1[i]));
                VOCPressureDic.Add(line1[0], newd);
                line = read.ReadLine();
            }
        }
        static double CalculatePearsonCorrelation(List<double> list1, List<double> list2)
        {
            // 确保两个列表长度相同  
            if (list1.Count != list2.Count)
            {
                throw new ArgumentException("Both lists must have the same length.");
            }

            int n = list1.Count;

            // 计算平均值  
            double mean1 = list1.Average();
            double mean2 = list2.Average();

            // 计算Pearson相关性系数  
            double numerator = 0; // 分子  
            double denominator1 = 0; // 分母的第一部分  
            double denominator2 = 0; // 分母的第二部分  

            for (int i = 0; i < n; i++)
            {
                double diff1 = list1[i] - mean1;
                double diff2 = list2[i] - mean2;

                numerator += diff1 * diff2;
                denominator1 += Math.Pow(diff1, 2);
                denominator2 += Math.Pow(diff2, 2);
            }

            // 防止除以0的情况  
            if (denominator1 == 0 || denominator2 == 0)
            {
                return 0; // 如果其中一个列表的所有值都相同，则相关性为0  
            }

            // 计算最终的Pearson相关性系数  
            double correlation = numerator / (Math.Sqrt(denominator1) * Math.Sqrt(denominator2));

            if (correlation < 0) correlation = 0;
            return correlation;
        }
        static void AntigenDistribution(string fold)
        {
            int i, j, k;
            StreamWriter writeMerged = new StreamWriter(fold + "/VOCWeeklyAntigenEpitopeT1.merged.txt");
            writeMerged.WriteLine("VOC\tDate\tWeeklyCases\tInfectionProp\tAntigenicDistanceT1\tNTD\tRBD\tSpikeOther\tEp01\tEp02\tEp03\tEp04\tEp05\tEp06\tEp07\tEp08\tEp09\tEp10\tEp11\tEp12\tRBDOther\tPressureCorrelation");
            for (i = 0; i < VOClist.Count(); i++)
            {
                for (j = 0; j < VOCDic[VOClist[i]].Weeklist.Count(); j++)
                {
                    string output = VOClist[i] + "\t" + VOCDic[VOClist[i]].Weeklist[j] + "\t" + VOCDic[VOClist[i]].WeekCases[j] + "\t" + VOCDic[VOClist[i]].WeekInfectionProp[j];
                    double weeklyAntigenicDistanceT1 = 0, weeklyNTD = 0, weeklyRBD = 0, weeklySpikeOther = 0, weeklyRBDOther = 0;
                    List<double> weeklyEpitope = new List<double>();
                    for (k = 0; k < 12; k++) weeklyEpitope.Add(0);

                    for (k = 0; k < VOCDic[VOClist[i]].SpikeMutList.Count(); k++)
                    {
                        if (VOCDic[VOClist[i]].SpikeMutList[k] == "")
                            continue;

                        string mutk = VOCDic[VOClist[i]].SpikeMutList[k].Substring(1, VOCDic[VOClist[i]].SpikeMutList[k].Length - 1);
                        int mutpos = Convert.ToInt32(mutk.Substring(0,mutk.Length-1));
                        string keyT1 = VOCDic[VOClist[i]].Weeklist[0] + "_" + VOCDic[VOClist[i]].SpikeMutList[k];//以第一个window为参考系
                        double freqT1 = VOCDic[VOClist[i]].DateMutFreqDic[keyT1];
                        string keyTn = VOCDic[VOClist[i]].Weeklist[j] + "_" + VOCDic[VOClist[i]].SpikeMutList[k];
                        double freqTn = VOCDic[VOClist[i]].DateMutFreqDic[keyTn];
                        string keyTn_1;
                        if (j > 0)
                            keyTn_1 = VOCDic[VOClist[i]].Weeklist[j - 1] + "_" + VOCDic[VOClist[i]].SpikeMutList[k];
                        else
                            keyTn_1 = VOCDic[VOClist[i]].Weeklist[0] + "_" + VOCDic[VOClist[i]].SpikeMutList[k];
                        double freqTn_1 = VOCDic[VOClist[i]].DateMutFreqDic[keyTn_1];

                        double distanceT1Freq = Math.Abs(freqTn - freqT1);
                        //double distanceT1Freq = Math.Abs(freqTn - freqTn_1);

                        //T1
                        if (AntigenicDistance.ContainsKey(mutk))
                        {
                            double dis = AntigenicDistance[mutk] * distanceT1Freq;
                            weeklyAntigenicDistanceT1 += dis;
                            if (mutpos >= 13 && mutpos <= 305)
                                weeklyNTD += dis;
                            else
                            if (mutpos >= 331 && mutpos <= 530)
                            {
                                double epiindex = 0;
                                weeklyRBD += dis;
                                if ("WT,Alpha,Delta".Contains(VOClist[i]))//preomicron
                                {
                                    for (int x = 0; x < PreOmicronEpitopeDic[mutpos].Count(); x++)
                                    {
                                        weeklyEpitope[x] += PreOmicronEpitopeDic[mutpos][x] * dis;
                                        epiindex += PreOmicronEpitopeDic[mutpos][x];
                                    }
                                    if (epiindex == 0)
                                        weeklyRBDOther += dis;
                                }
                                else//omicron
                                {
                                    for (int x = 0; x < OmicronEpitopeDic[mutpos].Count(); x++)
                                    {
                                        weeklyEpitope[x] += OmicronEpitopeDic[mutpos][x] * dis;
                                        epiindex += OmicronEpitopeDic[mutpos][x];
                                    }
                                    if (epiindex == 0)
                                        weeklyRBDOther += dis;
                                }
                            }
                            else
                                weeklySpikeOther += dis;
                        }
                    }
                    //"VOC\tDate\tWeeklyCases\tInfectionProp\tAntigenicDistanceT1\tNTD\tRBD\tSpikeOther\tEp01\tEp02\tEp03\tEp04\tEp05\tEp06\tEp07\tEp08\tEp09\tEp10\tEp11\tEp12\tRBDOther"
                    output += "\t" + Convert.ToString(weeklyAntigenicDistanceT1);
                    output += "\t" + Convert.ToString(weeklyNTD);
                    output += "\t" + Convert.ToString(weeklyRBD);
                    output += "\t" + Convert.ToString(weeklySpikeOther);
                    for(k=0;k<12;k++)
                        output += "\t" + Convert.ToString(weeklyEpitope[k]);            
                    output += "\t" + Convert.ToString(weeklyRBDOther);
                    output += "\t" + Convert.ToString(CalculatePearsonCorrelation(weeklyEpitope, VOCPressureDic[VOClist[i]]));
                    writeMerged.WriteLine(output);
                }
            }
            writeMerged.Close();
            return;
        }
        static void CalculateHaplotypeDistance(string fold)//计算单倍型到对应窗口其他序列的距离，为下一步计算相对增长优势做准备
        {
            int i, j, k;
            string[] filePaths = Directory.GetFiles(fold, "*.tsv", SearchOption.TopDirectoryOnly);
            for (i=0;i<VOClist.Count();i++)
            {
                for(j=0;j<filePaths.Count();j++)
                    if (filePaths[j].Contains(VOClist[i].Replace('/','_')+".Top"))
                    {
                        StreamReader read = new StreamReader(filePaths[j]);
                        StreamWriter write = new StreamWriter(filePaths[j] + ".Distance");
                        write.WriteLine("VOC\tHaploName\tDate\tDateNum\tHaploCount\tOtherCount\tHaploCases\tOtherCases\tAntigenicDistanceABS\tAntigenicDistance\tChangeACE2binding\tChangeCellEntry\tChangeImmuneEscape\tInfectionProportion");
                        string line = read.ReadLine();
                        line = read.ReadLine();
                        int lineind = 1;
                        while(line!=null)
                        {
                            double weeklyAntigenicDistanceABS = 0, weeklyAntigenicDistance = 0, weeklySpikeCellEntry = 0, weeklyRBDACE2binding = 0, weeklyEscapeScore = 0;
                            string[] line1 = line.Split('\t');
                            string[] mut1 = line1[3].Split(' ');
                            for (k = 0; k < VOCDic[VOClist[i]].SpikeMutList.Count();k++)
                            {
                                string mut = VOCDic[VOClist[i]].SpikeMutList[k];
                                string mutk = "";
                                if(mut!="") mutk = VOCDic[VOClist[i]].SpikeMutList[k].Substring(1, VOCDic[VOClist[i]].SpikeMutList[k].Length - 1);
                                double freq = VOCDic[VOClist[i]].DateMutFreqDic[line1[2] + "_" + mut];
                                if (mut1.Contains(mut))//计算差值
                                    freq = 1 - (freq - Convert.ToDouble(line1[5])) / (1 - Convert.ToDouble(line1[5]));//freq = 1 - freq;
                                else
                                    freq = 0 - freq;

                                if (AntigenicDistance.ContainsKey(mutk))
                                {
                                    weeklyAntigenicDistanceABS += AntigenicDistance[mutk] * Math.Abs(freq);
                                    weeklyAntigenicDistance += AntigenicDistance[mutk] * freq;
                                }
                                if (VOCDic[VOClist[i]].Spike_CellEntry.ContainsKey(mutk))
                                    weeklySpikeCellEntry += VOCDic[VOClist[i]].Spike_CellEntry[mutk] * freq;
                                if (VOCDic[VOClist[i]].RBD_ACE2Binding.ContainsKey(mutk))
                                    weeklyRBDACE2binding += VOCDic[VOClist[i]].RBD_ACE2Binding[mutk] * freq;
                                if (VOCDic[VOClist[i]].RBD_EscapeScore.ContainsKey(mutk))
                                    weeklyEscapeScore += VOCDic[VOClist[i]].RBD_EscapeScore[mutk] * freq;
                            }

                            string output = VOClist[i];
                            output += "\t" + line1[1];
                            output += "\t" + line1[2];
                            output += "\t" + Convert.ToString(lineind);
                            output += "\t" + line1[8];
                            output += "\t" + line1[9];
                            output += "\t" + line1[6];
                            output += "\t" + Convert.ToString(Convert.ToDouble(line1[7])- Convert.ToDouble(line1[6]));
                            output += "\t" + Convert.ToString(weeklyAntigenicDistanceABS);
                            output += "\t" + Convert.ToString(weeklyAntigenicDistance);
                            output += "\t" + Convert.ToString(weeklyRBDACE2binding);
                            output += "\t" + Convert.ToString(weeklySpikeCellEntry);
                            output += "\t" + Convert.ToString(weeklyEscapeScore);
                            output += "\t" + Convert.ToString(VOCDic[VOClist[i]].WeekInfectionProp[VOCDic[VOClist[i]].Weeklist.IndexOf(line1[2])]);
                            write.WriteLine(output);

                            line = read.ReadLine();
                            lineind++;
                        }
                        write.Close();
                    }
            }
            return;
        }
        static void CalculateHaplotypePreviousVOCDistance(string fold)//计算单倍型到前一个VOC对应窗口的距离，为下一步计算相对增长优势做准备
        {
            int i, j, k;
            Dictionary<string, string> prevocDic = new Dictionary<string, string>();
            prevocDic.Add("WT", "WT"); prevocDic.Add("Alpha", "WT"); prevocDic.Add("Delta", "Alpha"); prevocDic.Add("BA.1", "Delta");
            prevocDic.Add("BA.2", "BA.1"); prevocDic.Add("BA.4/5", "BA.2"); prevocDic.Add("BQ.1", "BA.4/5"); prevocDic.Add("XBB", "BQ.1");
            prevocDic.Add("BA.2.86", "XBB");
            string[] filePaths = Directory.GetFiles(fold, "*.tsv", SearchOption.TopDirectoryOnly);
            for (i = 0; i < VOClist.Count(); i++)
            {
                for (j = 0; j < filePaths.Count(); j++)
                    if (filePaths[j].Contains(VOClist[i].Replace('/', '_') + ".Top"))
                    {
                        StreamReader read = new StreamReader(filePaths[j]);
                        StreamWriter write = new StreamWriter(filePaths[j] + ".PreVOCDis");
                        write.WriteLine("VOC\tHaploName\tDate\tDateNum\tHaploCount\tOtherCount\tHaploCases\tOtherCases\tAntigenicDistanceABS\tAntigenicDistance\tChangeACE2binding\tChangeCellEntry\tChangeImmuneEscape\tInfectionProportion_PreVOC");
                        string line = read.ReadLine();
                        line = read.ReadLine();
                        int lineind = 1;
                        while (line != null)
                        {
                            double weeklyAntigenicDistanceABS = 0, weeklyAntigenicDistance = 0, weeklySpikeCellEntry = 0, weeklyRBDACE2binding = 0, weeklyEscapeScore = 0;
                            string[] line1 = line.Split('\t');
                            string[] mut1 = line1[3].Split(' ');
                            string prevoc = prevocDic[VOClist[i]];

                            if (VOCDic[prevoc].Weeklist.Contains(line1[2]))
                            {
                                for (k = 0; k < VOCDic[prevoc].SpikeMutList.Count(); k++)
                                {
                                    string mut = VOCDic[prevoc].SpikeMutList[k];
                                    string mutk = "";
                                    if (mut != "") mutk = VOCDic[prevoc].SpikeMutList[k].Substring(1, VOCDic[prevoc].SpikeMutList[k].Length - 1);
                                    double freq = VOCDic[prevoc].DateMutFreqDic[line1[2] + "_" + mut];
                                    if (mut1.Contains(mut))//计算差值
                                        freq = 1 - freq;//freq = 1 - freq;
                                    else
                                        freq = 0 - freq;

                                    if (AntigenicDistance.ContainsKey(mutk))
                                    {
                                        weeklyAntigenicDistanceABS += AntigenicDistance[mutk] * Math.Abs(freq);
                                        weeklyAntigenicDistance += AntigenicDistance[mutk] * freq;
                                    }
                                    if (VOCDic[prevoc].Spike_CellEntry.ContainsKey(mutk))
                                        weeklySpikeCellEntry += VOCDic[prevoc].Spike_CellEntry[mutk] * freq;
                                    if (VOCDic[prevoc].RBD_ACE2Binding.ContainsKey(mutk))
                                        weeklyRBDACE2binding += VOCDic[prevoc].RBD_ACE2Binding[mutk] * freq;
                                    if (VOCDic[prevoc].RBD_EscapeScore.ContainsKey(mutk))
                                        weeklyEscapeScore += VOCDic[prevoc].RBD_EscapeScore[mutk] * freq;
                                }

                                string output = prevoc;
                                output += "\t" + line1[1];
                                output += "\t" + line1[2];
                                output += "\t" + Convert.ToString(lineind);
                                output += "\t" + line1[8];
                                output += "\t" + line1[10];
                                output += "\t" + line1[6];
                                output += "\t" + Convert.ToString(Convert.ToDouble(line1[7]) - Convert.ToDouble(line1[6]));
                                output += "\t" + Convert.ToString(weeklyAntigenicDistanceABS);
                                output += "\t" + Convert.ToString(weeklyAntigenicDistance);
                                output += "\t" + Convert.ToString(weeklyRBDACE2binding);
                                output += "\t" + Convert.ToString(weeklySpikeCellEntry);
                                output += "\t" + Convert.ToString(weeklyEscapeScore);
                                output += "\t" + Convert.ToString(VOCDic[prevoc].WeekInfectionProp[VOCDic[prevoc].Weeklist.IndexOf(line1[2])]);
                                write.WriteLine(output);
                            }
                            line = read.ReadLine();
                            lineind++;
                        }
                        write.Close();
                    }
            }
            return;
        }
        static void CalculateHaplotypeT1Distance(string fold)//计算单倍型到对应VOC T1窗口的距离
        {
            int i, j, k;
            string[] filePaths = Directory.GetFiles(fold, "*.tsv", SearchOption.TopDirectoryOnly);
            for (i = 0; i < VOClist.Count(); i++)
            {
                for (j = 0; j < filePaths.Count(); j++)
                    if (filePaths[j].Contains(VOClist[i].Replace('/', '_') + ".Top"))
                    {
                        StreamReader read = new StreamReader(filePaths[j]);
                        StreamWriter write = new StreamWriter(filePaths[j] + ".T1Distance");
                        write.WriteLine("VOC\tHaploName\tDate\tDateNum\tHaploCount\tOtherCount\tHaploCases\tOtherCases\tAntigenicDistanceABS\tAntigenicDistance\tChangeACE2binding\tChangeCellEntry\tChangeImmuneEscape\tInfectionProportion");
                        string line = read.ReadLine();
                        line = read.ReadLine();
                        int lineind = 1;
                        string T1Date = "#TBD";
                        double weeklyAntigenicDistanceABS = 0, weeklyAntigenicDistance = 0, weeklySpikeCellEntry = 0, weeklyRBDACE2binding = 0, weeklyEscapeScore = 0;
                        while (line != null)
                        {
                            
                            string[] line1 = line.Split('\t');
                            string[] mut1 = line1[3].Split(' ');
                            if (T1Date == "#TBD")
                            {
                                T1Date = line1[2];
                                for (k = 0; k < VOCDic[VOClist[i]].SpikeMutList.Count(); k++)
                                {
                                    string mut = VOCDic[VOClist[i]].SpikeMutList[k];
                                    string mutk = "";
                                    if (mut != "") mutk = VOCDic[VOClist[i]].SpikeMutList[k].Substring(1, VOCDic[VOClist[i]].SpikeMutList[k].Length - 1);
                                    double freq = VOCDic[VOClist[i]].DateMutFreqDic[T1Date + "_" + mut];
                                    if (mut1.Contains(mut))//计算差值
                                        freq = 1 - freq;//freq = 1 - freq;
                                    else
                                        freq = 0 - freq;

                                    if (AntigenicDistance.ContainsKey(mutk))
                                    {
                                        weeklyAntigenicDistanceABS += AntigenicDistance[mutk] * Math.Abs(freq);
                                        weeklyAntigenicDistance += AntigenicDistance[mutk] * freq;
                                    }
                                    if (VOCDic[VOClist[i]].Spike_CellEntry.ContainsKey(mutk))
                                        weeklySpikeCellEntry += VOCDic[VOClist[i]].Spike_CellEntry[mutk] * freq;
                                    if (VOCDic[VOClist[i]].RBD_ACE2Binding.ContainsKey(mutk))
                                        weeklyRBDACE2binding += VOCDic[VOClist[i]].RBD_ACE2Binding[mutk] * freq;
                                    if (VOCDic[VOClist[i]].RBD_EscapeScore.ContainsKey(mutk))
                                        weeklyEscapeScore += VOCDic[VOClist[i]].RBD_EscapeScore[mutk] * freq;
                                }
                            }

                            string output = VOClist[i];
                            output += "\t" + line1[1];
                            output += "\t" + line1[2];
                            output += "\t" + Convert.ToString(lineind);
                            output += "\t" + line1[8];
                            output += "\t" + line1[9];
                            output += "\t" + line1[6];
                            output += "\t" + Convert.ToString(Convert.ToDouble(line1[7]) - Convert.ToDouble(line1[6]));
                            output += "\t" + Convert.ToString(weeklyAntigenicDistanceABS);
                            output += "\t" + Convert.ToString(weeklyAntigenicDistance);
                            output += "\t" + Convert.ToString(weeklyRBDACE2binding);
                            output += "\t" + Convert.ToString(weeklySpikeCellEntry);
                            output += "\t" + Convert.ToString(weeklyEscapeScore);
                            output += "\t" + Convert.ToString(VOCDic[VOClist[i]].WeekInfectionProp[VOCDic[VOClist[i]].Weeklist.IndexOf(line1[2])]);
                            write.WriteLine(output);

                            line = read.ReadLine();
                            lineind++;
                        }
                        write.Close();
                    }
            }
            return;
        }
        static void CalculateHaplotypePreviousVOCT1Distance(string fold)//计算单倍型到前一个VOCT1窗口的距离
        {
            int i, j, k;
            Dictionary<string, string> prevocDic = new Dictionary<string, string>();
            prevocDic.Add("WT", "WT"); prevocDic.Add("Alpha", "WT"); prevocDic.Add("Delta", "Alpha"); prevocDic.Add("BA.1", "Delta");
            prevocDic.Add("BA.2", "BA.1"); prevocDic.Add("BA.4/5", "BA.2"); prevocDic.Add("BQ.1", "BA.4/5"); prevocDic.Add("XBB", "BQ.1");
            prevocDic.Add("BA.2.86", "XBB");
            string[] filePaths = Directory.GetFiles(fold, "*.Top1.tsv", SearchOption.TopDirectoryOnly);
            //VOClist = "WT,Alpha,Delta,BA.1,BA.2,BA.4/5,BQ.1,XBB,BA.2.86".Split(',');
            VOClist = "BA.2".Split(',');
            for (i = 0; i < VOClist.Count(); i++)
            {
                for (j = 0; j < filePaths.Count(); j++)
                    if (filePaths[j].Contains(VOClist[i].Replace('/', '_') + ".Top"))
                    {
                        StreamReader read = new StreamReader(filePaths[j]);
                        StreamWriter write = new StreamWriter(filePaths[j] + ".PreVOCT1Dis");
                        write.WriteLine("VOC\tHaploName\tDate\tDateNum\tHaploCount\tOtherCount\tHaploCases\tOtherCases\tAntigenicDistanceABS\tAntigenicDistance\tChangeACE2binding\tChangeCellEntry\tChangeImmuneEscape\tInfectionProportion_PreVOC\tInfectionProportion");
                        string line = read.ReadLine();
                        line = read.ReadLine();
                        int lineind = 1;
                        while (line != null)
                        {
                            double weeklyAntigenicDistanceABS = 0, weeklyAntigenicDistance = 0, weeklySpikeCellEntry = 0, weeklyRBDACE2binding = 0, weeklyEscapeScore = 0;
                            string[] line1 = line.Split('\t');
                            string[] mut1 = line1[3].Split(' ');
                            string prevoc = prevocDic[VOClist[i]];

                            if (true)
                            {
                                for (k = 0; k < VOCDic[prevoc].SpikeMutList.Count(); k++)
                                {
                                    string mut = VOCDic[prevoc].SpikeMutList[k];
                                    string mutk = "";
                                    if (mut != "") mutk = VOCDic[prevoc].SpikeMutList[k].Substring(1, VOCDic[prevoc].SpikeMutList[k].Length - 1);
                                    double freq = VOCDic[prevoc].DateMutFreqDic[VOCDic[prevoc].Weeklist[0] + "_" + mut];
                                    if (mut1.Contains(mut))//计算差值
                                        freq = 1 - freq;//freq = 1 - freq;
                                    else
                                        freq = 0 - freq;

                                    if (AntigenicDistance.ContainsKey(mutk))
                                    {
                                        weeklyAntigenicDistanceABS += AntigenicDistance[mutk] * Math.Abs(freq);
                                        weeklyAntigenicDistance += AntigenicDistance[mutk] * freq;
                                    }
                                    if (VOCDic[prevoc].Spike_CellEntry.ContainsKey(mutk))
                                        weeklySpikeCellEntry += VOCDic[prevoc].Spike_CellEntry[mutk] * freq;
                                    if (VOCDic[prevoc].RBD_ACE2Binding.ContainsKey(mutk))
                                        weeklyRBDACE2binding += VOCDic[prevoc].RBD_ACE2Binding[mutk] * freq;
                                    if (VOCDic[prevoc].RBD_EscapeScore.ContainsKey(mutk))
                                        weeklyEscapeScore += VOCDic[prevoc].RBD_EscapeScore[mutk] * freq;
                                }

                                string output = prevoc;
                                output += "\t" + line1[1];
                                output += "\t" + line1[2];
                                output += "\t" + Convert.ToString(lineind);
                                output += "\t" + line1[8];
                                output += "\t" + line1[10];
                                output += "\t" + line1[6];
                                output += "\t" + Convert.ToString(Convert.ToDouble(line1[7]) - Convert.ToDouble(line1[6]));
                                output += "\t" + Convert.ToString(weeklyAntigenicDistanceABS);
                                output += "\t" + Convert.ToString(weeklyAntigenicDistance);
                                output += "\t" + Convert.ToString(weeklyRBDACE2binding);
                                output += "\t" + Convert.ToString(weeklySpikeCellEntry);
                                output += "\t" + Convert.ToString(weeklyEscapeScore);
                                output += "\tNA";// + Convert.ToString(VOCDic[prevoc].WeekInfectionProp[VOCDic[prevoc].Weeklist.IndexOf(line1[2])]);
                                output += "\tNA";// + Convert.ToString(VOCDic[VOClist[i]].WeekInfectionProp[VOCDic[VOClist[i]].Weeklist.IndexOf(line1[2])]);
                                write.WriteLine(output);
                            }
                            line = read.ReadLine();
                            lineind++;
                        }
                        write.Close();
                    }
            }
            return;
        }
        static void Main(string[] args)
        {
            ReadinDMS();//read processed DMS data
            ReadinMatrix();//read weekly mutation frequency matrix
            CalculateDistance("//NAS8500/g/VariationMutation/AiDMS/Spike0417/WeeklyMatrix");//calculate weekly distance
            //CalculateAntigenDistanceMatrix("//NAS8500/g/VariationMutation/AiDMS/Spike0417/WeeklyMatrix/MDS");//calcualte pairwise week distance
            //ReadinEpitope();//read RBD epitope information
            //AntigenDistribution("//NAS8500/g/VariationMutation/AiDMS/Spike0417/WeeklyMatrix/Epitope");
            CalculateHaplotypeDistance("//NAS8500/g/VariationMutation/AiDMS/Spike0417/Diversity/WeeklyMax/Haplotypes");//Calculate the distance from haplotypes to the remaining sequences in the corresponding window, in preparation for the next step of calculating relative growth advantages
            //CalculateHaplotypeT1Distance("//NAS8500/g/VariationMutation/AiDMS/Spike0417/Diversity/WeeklyMax/Haplotypes");//
            CalculateHaplotypePreviousVOCDistance("//NAS8500/g/VariationMutation/AiDMS/Spike0417/Diversity/HaploCasesTop1");//Calculate the distance from the haplotype to the overlapping window of the previous VOC, in preparation for the next step of calculating relative growth advantages
            //CalculateHaplotypePreviousVOCT1Distance("//NAS8500/g/VariationMutation/AiDMS/Spike0417/Diversity/WeeklyMax/Haplotypes");//
            return;
        }
    }
}
