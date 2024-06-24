using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Xml.Schema;

namespace 新_计算每周的haplo分布关联感染人数
{
    public class SpikeHaplotype
    {
        public Dictionary<string, char> SpikeMutDic = new Dictionary<string, char>();//amino acid mutations in spike
        public Dictionary<string, int> LineageCount = new Dictionary<string, int>();
        public int totalSeq = 0;
        public string representingLineage = "TBD";//choose the one with most sequence

    }
    public class Sequence
    {
        public string Seqid;//accesion id
        public string SpikeHaploMut; //use mutation as dictionary key
        public string CollectionDate;
        public string CollectionLocation;
        public string Lineage;
        public string VOCGroup = "Other";
    }
    public class Week
    {
        public string Weekdate; //DateTime - 7
        public int Weekstart;
        public int Weekend;
        public double Cases;
        public double CasesSmoothed;
        public List<Sequence> WeekSeqList = new List<Sequence>();

        public Dictionary<string, int> WeekVOCCountDic = new Dictionary<string, int>();
        public Dictionary<string, Dictionary<string, int>> WeekVOCHaploCountDic = new Dictionary<string, Dictionary<string, int>>();
        
        public Dictionary<string, double> WeekVOCGenomeNUC = new Dictionary<string, double>();
        public Dictionary<string, double> WeekVOCGenomeCount = new Dictionary<string, double>();
    }
    internal class Program
    {
        static List<string> VOCList = new List<string>();
        static Dictionary<string, int> VOCStartDate = new Dictionary<string, int>();//date start analysis
        static Dictionary<string, int> VOCEndDate = new Dictionary<string, int>();//date end analysis
        static Dictionary<string, double> VOCTotalCases = new Dictionary<string, double>();

        static Dictionary<string, string> Lineage2VOC = new Dictionary<string, string>();//lineage到VOC的对照表
        static List<Week> Weeklist = new List<Week>();//存储所有的周
        static Dictionary<string, SpikeHaplotype> HaploDic = new Dictionary<string, SpikeHaplotype>();//从突变到haplo的对照表

        static List<string> SpikeMutList = new List<string>();//记录所有的Spike AA突变 输出矩阵

        static void ReadinVOCTimeWindow(string file)//读入lineage流行的起止时间
        {
            int i, j, k;
            StreamReader read = new StreamReader(file);
            string line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                VOCStartDate.Add(line1[0], Convert.ToInt32(line1[1]));
                VOCEndDate.Add(line1[0], Convert.ToInt32(line1[2]));
                line = read.ReadLine();
            }
            read.Close();
            return;
        }
        static void ReadingLineagenote(string file)//读入lineage对照表
        {
            StreamReader read = new StreamReader(file);
            string line = read.ReadLine();
            line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                if (line1[2] != "Other")
                {
                    if (!VOCList.Contains(line1[2]))
                        VOCList.Add(line1[2]);
                    Lineage2VOC.Add(line1[0], line1[2]);
                }
                line = read.ReadLine();
            }
            read.Close();
            return;
        }
        static void ReadinCases(string file, string country)//读入每周感染人数，构建time window
        {
            int i, j, k;
            int datestart = -1;
            StreamReader read = new StreamReader(file);
            string line = read.ReadLine();
            line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split(',');
                if (line1[2]==country)
                {
                    if(datestart!=-1)
                    {
                        Week neww = new Week();
                        neww.Weekstart = datestart;
                        neww.Weekend = Convert.ToInt32(line1[0].Substring(0, 4)) * 10000 + Convert.ToInt32(line1[0].Substring(5, 2)) * 100 + Convert.ToInt32(line1[0].Substring(8, 2));
                        neww.Weekdate = line1[0];
                        for (i = 0; i < VOCList.Count(); i++)
                        {
                            neww.WeekVOCCountDic.Add(VOCList[i], 0);
                            neww.WeekVOCGenomeNUC.Add(VOCList[i], 0);
                            neww.WeekVOCGenomeCount.Add(VOCList[i], 0);
                        }

                        if (line1[4] == "") line1[4] = "0";
                        neww.Cases = Convert.ToDouble(line1[4]);
                        neww.CasesSmoothed = Convert.ToDouble(line1[4]);
                        Weeklist.Add(neww);
                    }
                    datestart = Convert.ToInt32(line1[0].Substring(0, 4)) * 10000 + Convert.ToInt32(line1[0].Substring(5, 2)) * 100 + Convert.ToInt32(line1[0].Substring(8, 2));
                }
                line = read.ReadLine();
            }
            read.Close();

            //滑动平均
            for (i = 1; i < Weeklist.Count() - 1; i++)
                Weeklist[i].CasesSmoothed = (Weeklist[i-1].Cases + Weeklist[i].Cases + Weeklist[i+1].Cases) / 3;
            return;
        }
        static void ReadinSeq(string fold, string country)//读入序列，按采样时间扔到time window里
        {
            int i, j, k;
            StreamReader read = new StreamReader(fold + "/sample_mut_loc_time." + country + ".tsv.SpikeAAmut");
            string line = read.ReadLine();
            while(line!=null)
            {
                string[] line1 = line.Split('\t');
                if (line1[3].Length == 10)
                {
                    Sequence news = new Sequence();
                    news.Seqid = line1[0];
                    news.CollectionLocation = line1[2];
                    news.CollectionDate = line1[3];
                    news.Lineage = line1[4];
                    news.SpikeHaploMut = line1[1];
                    if (Lineage2VOC.ContainsKey(line1[4]))
                        news.VOCGroup = Lineage2VOC[line1[4]];

                    if (HaploDic.ContainsKey(line1[1]))
                    {
                        if (!HaploDic[line1[1]].LineageCount.ContainsKey(line1[4]))
                            HaploDic[line1[1]].LineageCount.Add(line1[4], 0);
                        HaploDic[line1[1]].LineageCount[line1[4]]++;
                    }
                    else
                    {
                        SpikeHaplotype newhap = new SpikeHaplotype();
                        newhap.LineageCount.Add(line1[4], 1);
                        string[] mut1 = line1[1].Split(' ');
                        for (i = 0; i < mut1.Count(); i++)
                            newhap.SpikeMutDic.Add(mut1[i], '1');
                        HaploDic.Add(line1[1], newhap);
                    }
                    HaploDic[line1[1]].totalSeq++;

                    int date = Convert.ToInt32(line1[3].Substring(0, 4)) * 10000 + Convert.ToInt32(line1[3].Substring(5, 2)) * 100 + Convert.ToInt32(line1[3].Substring(8, 2));
                    for (i = 0; i < Weeklist.Count(); i++)
                        if (Weeklist[i].Weekstart <= date && date < Weeklist[i].Weekend)
                        {
                            Weeklist[i].WeekSeqList.Add(news);
                            if (news.VOCGroup != "Other")
                                Weeklist[i].WeekVOCCountDic[news.VOCGroup]++;
                        }
                }
                line = read.ReadLine();
            }
            read.Close();

            read = new StreamReader(fold + "/sample_mut_loc_time.tsv");
            line = read.ReadLine();
            while (line != null)
            {
                string[] line1 = line.Split('\t');
                if (line1[3].Length == 10)
                {
                    int date = Convert.ToInt32(line1[3].Substring(0, 4)) * 10000 + Convert.ToInt32(line1[3].Substring(5, 2)) * 100 + Convert.ToInt32(line1[3].Substring(8, 2));
                    for (i = 0; i < Weeklist.Count(); i++)
                        if (Weeklist[i].Weekstart <= date && date < Weeklist[i].Weekend)
                        {
                            if (Lineage2VOC.ContainsKey(line1[4]) && Lineage2VOC[line1[4]] != "Other")
                            {
                                string[] mut1 = line1[1].Split(' ');
                                Weeklist[i].WeekVOCGenomeNUC[Lineage2VOC[line1[4]]] += mut1.Count();
                                Weeklist[i].WeekVOCGenomeCount[Lineage2VOC[line1[4]]]++;
                            }
                        }
                }
                line = read.ReadLine();
            }
            read.Close();

            return;
        }
        static void CalculateRepresentingLineage()//计算每个haplo的代表性lineage
        {
            int i, j, k;
            foreach(string val in HaploDic.Keys)
            {
                int maxcount = 0;
                foreach(string vallin in HaploDic[val].LineageCount.Keys)
                    if (HaploDic[val].LineageCount[vallin] > maxcount)
                    {
                        maxcount = HaploDic[val].LineageCount[vallin];
                        HaploDic[val].representingLineage = vallin;
                    }
            }
            return;
        }
        static void CalculateVOCCases(string fold, string country)//根据序列和感染人数计算每个voc的感染人数
        {
            StreamWriter write = new StreamWriter(fold + "/VOCWeeklyCases." + country + ".tsv");
            StreamWriter writeSeqCount = new StreamWriter(fold + "/VOCWeeklySequence." + country + ".tsv");
            StreamWriter writePlot = new StreamWriter(fold + "/VOCWeeklyCases_ForPlot." + country + ".tsv");
            int i, j, k;
            string output = "Date";
            for (i = 0; i < VOCList.Count; i++)
            {
                output += "\t" + VOCList[i];
                VOCTotalCases.Add(VOCList[i], 0);
            }
            write.WriteLine(output);
            writeSeqCount.WriteLine(output);
            writePlot.WriteLine("Date\tVOC\tWeeklyCases");

            for (i=0;i<Weeklist.Count;i++)
            {
                output = Weeklist[i].Weekdate;
                for (j = 0; j < VOCList.Count; j++)
                {
                    VOCTotalCases[VOCList[j]] += Convert.ToDouble(Weeklist[i].WeekVOCCountDic[VOCList[j]]) / Weeklist[i].WeekSeqList.Count() * Weeklist[i].CasesSmoothed;
                    writePlot.WriteLine(Weeklist[i].Weekdate + "\t" + VOCList[j] +"\t" + Convert.ToString(Convert.ToDouble(Weeklist[i].WeekVOCCountDic[VOCList[j]]) / Weeklist[i].WeekSeqList.Count() * Weeklist[i].CasesSmoothed));
                    output += "\t" + Convert.ToString(Convert.ToDouble(Weeklist[i].WeekVOCCountDic[VOCList[j]]) / Weeklist[i].WeekSeqList.Count() * Weeklist[i].CasesSmoothed);
                }
                write.WriteLine(output);
                writePlot.WriteLine(Weeklist[i].Weekdate + "\tTotal\t" + Convert.ToString(Weeklist[i].CasesSmoothed));

                output = Weeklist[i].Weekdate;
                for (j = 0; j < VOCList.Count; j++)
                {
                    output += "\t" + Convert.ToString(Weeklist[i].WeekVOCCountDic[VOCList[j]]);
                }
                writeSeqCount.WriteLine(output); 
            }
            write.Close();
            writePlot.Close();
            writeSeqCount.Close();
            return;
        }
        static void CalculateHaploDiversity(string fold, string country)//计算单倍型diversity随时间的变化
        {
            int i, j, k;
            StreamWriter write = new StreamWriter(fold + "/VOCWeeklyHaplo." + country + ".tsv");
            write.WriteLine("VOC\tDate\tCases\tHaploNum\tHaploNumPer100Seq");
            for(i=0;i<VOCList.Count();i++)
            {
                
                
                for(j=0;j<Weeklist.Count;j++)
                {
                    if (Weeklist[j].Weekstart >= VOCStartDate[VOCList[i]] && Weeklist[j].Weekend <= VOCEndDate[VOCList[i]])
                    {
                        string output = VOCList[i] + "\t" + Weeklist[j].Weekdate;
                        output += "\t" + Convert.ToString(Convert.ToDouble(Weeklist[j].WeekVOCCountDic[VOCList[i]]) / Weeklist[j].WeekSeqList.Count() * Weeklist[j].CasesSmoothed);
                        Dictionary<string, int> tmpd = new Dictionary<string, int>();
                        for (k = 0; k < Weeklist[j].WeekSeqList.Count(); k++)
                            if (Weeklist[j].WeekSeqList[k].VOCGroup == VOCList[i])
                                if (!tmpd.ContainsKey(Weeklist[j].WeekSeqList[k].SpikeHaploMut))
                                    tmpd.Add(Weeklist[j].WeekSeqList[k].SpikeHaploMut, 1);
                        output += "\t" + Convert.ToString(tmpd.Count());
                        output += "\t" + Convert.ToString(Convert.ToDouble( tmpd.Count()) / Weeklist[j].WeekVOCCountDic[VOCList[i]]);
                        write.WriteLine(output);
                    }
                }
                
            }
            write.Close();
            return;
        }
        static void CalculateVOCHaplotype(string fold,string country)//计算每个VOC下haplotype的分布
        {
            StreamWriter write = new StreamWriter(fold + "/VOCHaploCases." + country + ".tsv");
            StreamWriter writeProp = new StreamWriter(fold + "/VOCHaploPropTopProp." + country + ".tsv");
            StreamWriter writeCase = new StreamWriter(fold + "/VOCHaploCaseTopProp." + country + ".tsv");
            write.WriteLine("VOC\tTop\tMut\tLineage\tCasesProportion");
            writeProp.WriteLine("VOC\tWeek\tMut\tLineage\tProportion");
            writeCase.WriteLine("VOC\tWeek\tMut\tLineage\tProportion");
            int i, j, k;
            for(i=0;i<Weeklist.Count();i++)
            {
                for(j=0;j<VOCList.Count();j++)
                {
                    Weeklist[i].WeekVOCHaploCountDic.Add(VOCList[j], new Dictionary<string, int>());
                    for (k = 0; k < Weeklist[i].WeekSeqList.Count();k++)
                    {
                        if (Weeklist[i].WeekSeqList[k].VOCGroup == VOCList[j])
                        {
                            if (Weeklist[i].WeekVOCHaploCountDic[VOCList[j]].ContainsKey(Weeklist[i].WeekSeqList[k].SpikeHaploMut))
                                Weeklist[i].WeekVOCHaploCountDic[VOCList[j]][Weeklist[i].WeekSeqList[k].SpikeHaploMut]++;
                            else
                                Weeklist[i].WeekVOCHaploCountDic[VOCList[j]].Add(Weeklist[i].WeekSeqList[k].SpikeHaploMut, 1);
                        }
                    }
                }
            }

            for (j = 0; j < VOCList.Count(); j++)
            {
                Dictionary<string, double> haploCaseDic = new Dictionary<string, double>();
                Dictionary<string, double> haploMaxPropDic = new Dictionary<string, double>();
                double totalVOCcases = 0;
                for (i = 0; i < Weeklist.Count(); i++)
                {
                    foreach (string val in Weeklist[i].WeekVOCHaploCountDic[VOCList[j]].Keys)
                    {
                        if (haploCaseDic.ContainsKey(val))//Cases
                        {
                            haploCaseDic[val] += Convert.ToDouble(Weeklist[i].WeekVOCHaploCountDic[VOCList[j]][val]) / Weeklist[i].WeekSeqList.Count() * Weeklist[i].CasesSmoothed;
                        }
                        else
                        {
                            haploCaseDic.Add(val, Convert.ToDouble(Weeklist[i].WeekVOCHaploCountDic[VOCList[j]][val]) / Weeklist[i].WeekSeqList.Count() * Weeklist[i].CasesSmoothed);
                        }

                        if (Weeklist[i].Weekstart >= VOCStartDate[VOCList[j]] && Weeklist[i].Weekend <= VOCEndDate[VOCList[j]])
                        {
                            if (haploMaxPropDic.ContainsKey(val))//Proportion
                            {
                                    if (Convert.ToDouble(Weeklist[i].WeekVOCHaploCountDic[VOCList[j]][val]) / Weeklist[i].WeekVOCCountDic[VOCList[j]] > haploMaxPropDic[val])
                                        haploMaxPropDic[val] = Convert.ToDouble(Weeklist[i].WeekVOCHaploCountDic[VOCList[j]][val]) / Weeklist[i].WeekVOCCountDic[VOCList[j]];
                            }
                            else
                            {
                                    haploMaxPropDic.Add(val, Convert.ToDouble(Weeklist[i].WeekVOCHaploCountDic[VOCList[j]][val]) / Weeklist[i].WeekVOCCountDic[VOCList[j]]);
                            }
                        }
                    }
                    totalVOCcases += Convert.ToDouble(Weeklist[i].WeekVOCCountDic[VOCList[j]]) / Weeklist[i].WeekSeqList.Count() * Weeklist[i].CasesSmoothed;
                }
                List<KeyValuePair<string, double>> sortedDictCases = new List<KeyValuePair<string, double>>(haploCaseDic.OrderBy(kvp => kvp.Value));
                List<KeyValuePair<string, double>> sortedDictSeqs = new List<KeyValuePair<string, double>>(haploMaxPropDic.OrderBy(kvp => kvp.Value));

                Dictionary<string, string> prevocDic = new Dictionary<string, string>();
                prevocDic.Add("WT", "WT"); prevocDic.Add("Alpha", "WT"); prevocDic.Add("Delta", "Alpha"); prevocDic.Add("BA.1", "Delta");
                prevocDic.Add("BA.2", "BA.1"); prevocDic.Add("BA.4/5", "BA.2"); prevocDic.Add("BQ.1", "BA.4/5"); prevocDic.Add("XBB", "BQ.1");
                prevocDic.Add("BA.2.86", "XBB");
                
                //输出感染人数前20的
                for (k = 1; k <= 2; k++)
                {
                    string output = VOCList[j] + "\t" + Convert.ToString(k);

                    if (k <= 20)
                    {
                        output += "\t" + Convert.ToString(sortedDictCases[sortedDictCases.Count - k].Key);
                        output += "\t" + Convert.ToString(HaploDic[sortedDictCases[sortedDictCases.Count - k].Key].representingLineage);
                        output += "\t" + Convert.ToString(sortedDictCases[sortedDictCases.Count - k].Value / totalVOCcases);
                        write.WriteLine(output);
                    }

                    if(k<=10)
                    for (i = 0; i < Weeklist.Count(); i++)
                    {
                        if (Weeklist[i].Weekstart >= VOCStartDate[VOCList[j]] && Weeklist[i].Weekend <= VOCEndDate[VOCList[j]])
                        {
                            output = VOCList[j] + "\t" + Weeklist[i].Weekdate + "\t" + sortedDictCases[sortedDictCases.Count - k].Key;
                            output += "\tTop" + Convert.ToString(k-1) + HaploDic[sortedDictCases[sortedDictCases.Count - k].Key].representingLineage;

                            if (Weeklist[i].WeekVOCHaploCountDic[VOCList[j]].ContainsKey(sortedDictCases[sortedDictCases.Count - k].Key))
                                output += "\t" + Convert.ToString(Convert.ToDouble(Weeklist[i].WeekVOCHaploCountDic[VOCList[j]][sortedDictCases[sortedDictCases.Count - k].Key]) / Weeklist[i].WeekVOCCountDic[VOCList[j]]);
                            else
                                output += "\t0";
                            writeCase.WriteLine(output);
                        }
                    }

                    if (k <= 80)
                    {
                        StreamWriter writehaplocase = new StreamWriter(fold + "/HaploCasesTop1/" + VOCList[j].Replace('/','_') + ".Top" + Convert.ToString(k) + ".tsv");
                        writehaplocase.WriteLine("VOC\tTopK\tWeek\tMut\tLineage\tProportion\tHaploWeeklyCases\tVOCWeeklyCases\tHaploCount\tOtherCount\tPreviousVOCCount");
                        for (i = 0; i < Weeklist.Count(); i++)
                        {
                            if (Weeklist[i].Weekstart >= VOCStartDate[VOCList[j]] && Weeklist[i].Weekend <= VOCEndDate[VOCList[j]])
                            {
                                output = VOCList[j] + "\t" + Convert.ToString(k) + "\t" + Weeklist[i].Weekdate + "\t" + sortedDictCases[sortedDictCases.Count - k].Key;
                                output += "\t" + HaploDic[sortedDictCases[sortedDictCases.Count - k].Key].representingLineage;

                                if (Weeklist[i].WeekVOCHaploCountDic[VOCList[j]].ContainsKey(sortedDictCases[sortedDictCases.Count - k].Key))
                                {
                                    output += "\t" + Convert.ToString(Convert.ToDouble(Weeklist[i].WeekVOCHaploCountDic[VOCList[j]][sortedDictCases[sortedDictCases.Count - k].Key]) / Weeklist[i].WeekVOCCountDic[VOCList[j]]);//haplo proportion
                                    output += "\t" + Convert.ToString(Convert.ToDouble(Weeklist[i].WeekVOCHaploCountDic[VOCList[j]][sortedDictCases[sortedDictCases.Count - k].Key]) / Weeklist[i].WeekSeqList.Count() * Weeklist[i].CasesSmoothed);//haplo cases
                                    output += "\t" + Convert.ToString(Convert.ToDouble(Weeklist[i].WeekVOCCountDic[VOCList[j]]) / Weeklist[i].WeekSeqList.Count() * Weeklist[i].CasesSmoothed);//VOC cases
                                    output += "\t" + Convert.ToString(Weeklist[i].WeekVOCHaploCountDic[VOCList[j]][sortedDictCases[sortedDictCases.Count - k].Key]);//haplo count
                                    output += "\t" + Convert.ToString(Weeklist[i].WeekVOCCountDic[VOCList[j]] - Weeklist[i].WeekVOCHaploCountDic[VOCList[j]][sortedDictCases[sortedDictCases.Count - k].Key]);//other count
                                }
                                else
                                {
                                    output += "\t0\t0";
                                    output += "\t" + Convert.ToString(Convert.ToDouble(Weeklist[i].WeekVOCCountDic[VOCList[j]]) / Weeklist[i].WeekSeqList.Count() * Weeklist[i].CasesSmoothed);//voc cases
                                    output += "\t0";//haplo count
                                    output += "\t" + Convert.ToString(Weeklist[i].WeekVOCCountDic[VOCList[j]]);//other count
                                }
                                output += "\t" + Convert.ToString(Weeklist[i].WeekVOCCountDic[prevocDic[VOCList[j]]]);
                                writehaplocase.WriteLine(output);
                            }
                        }
                        writehaplocase.Close();
                    }
                }

                
                //输出最大占比
                for (k = 1; k < sortedDictSeqs.Count; k++)
                {
                    //if (haploMaxPropDic[sortedDictSeqs[sortedDictSeqs.Count - k].Key] >= 0.001 && HaploDic[sortedDictSeqs[sortedDictSeqs.Count - k].Key].totalSeq >= 10)
                    if (HaploDic[sortedDictSeqs[sortedDictSeqs.Count - k].Key].totalSeq >= 10)
                    {
                        
                        StreamWriter writehaploseqs = new StreamWriter(fold + "/WeeklyMax/Haplotypes/" + VOCList[j].Replace('/', '_') + ".Top" + Convert.ToString(k) + ".tsv");
                        writehaploseqs.WriteLine("VOC\tTopK\tWeek\tMut\tLineage\tProportion\tHaploWeeklyCases\tVOCWeeklyCases\tHaploCount\tOtherCount\tPreviousVOCCount");
                        for (i = 0; i < Weeklist.Count(); i++)
                        {
                            if (Weeklist[i].Weekstart >= VOCStartDate[VOCList[j]] && Weeklist[i].Weekend <= VOCEndDate[VOCList[j]])
                            {
                                string output = VOCList[j] + "\t" + Convert.ToString(k) + "\t" + Weeklist[i].Weekdate + "\t" + sortedDictSeqs[sortedDictSeqs.Count - k].Key;
                                output += "\t" + HaploDic[sortedDictSeqs[sortedDictSeqs.Count - k].Key].representingLineage;

                                if (Weeklist[i].WeekVOCHaploCountDic[VOCList[j]].ContainsKey(sortedDictSeqs[sortedDictSeqs.Count - k].Key))
                                {
                                    output += "\t" + Convert.ToString(Convert.ToDouble(Weeklist[i].WeekVOCHaploCountDic[VOCList[j]][sortedDictSeqs[sortedDictSeqs.Count - k].Key]) / Weeklist[i].WeekVOCCountDic[VOCList[j]]);//haplo proportion
                                    output += "\t" + Convert.ToString(Convert.ToDouble(Weeklist[i].WeekVOCHaploCountDic[VOCList[j]][sortedDictSeqs[sortedDictSeqs.Count - k].Key]) / Weeklist[i].WeekSeqList.Count() * Weeklist[i].CasesSmoothed);//haplo cases
                                    output += "\t" + Convert.ToString(Convert.ToDouble(Weeklist[i].WeekVOCCountDic[VOCList[j]]) / Weeklist[i].WeekSeqList.Count() * Weeklist[i].CasesSmoothed);//VOC cases
                                    output += "\t" + Convert.ToString(Weeklist[i].WeekVOCHaploCountDic[VOCList[j]][sortedDictSeqs[sortedDictSeqs.Count - k].Key]);//haplo count
                                    output += "\t" + Convert.ToString(Weeklist[i].WeekVOCCountDic[VOCList[j]] - Weeklist[i].WeekVOCHaploCountDic[VOCList[j]][sortedDictSeqs[sortedDictSeqs.Count - k].Key]);//other count
                                }
                                else
                                {
                                    output += "\t0\t0";
                                    output += "\t" + Convert.ToString(Convert.ToDouble(Weeklist[i].WeekVOCCountDic[VOCList[j]]) / Weeklist[i].WeekSeqList.Count() * Weeklist[i].CasesSmoothed);//voc cases
                                    output += "\t0";//haplo count
                                    output += "\t" + Convert.ToString(Weeklist[i].WeekVOCCountDic[VOCList[j]]);//other count
                                }
                                output += "\t" + Convert.ToString(Weeklist[i].WeekVOCCountDic[prevocDic[VOCList[j]]]);
                                writehaploseqs.WriteLine(output);
                            }
                        }
                        writehaploseqs.Close();
                    }
                }

                //输出最高频率大于【】的
                for (i = 0; i < Weeklist.Count(); i++)
                {
                    if (Weeklist[i].Weekstart >= VOCStartDate[VOCList[j]] && Weeklist[i].Weekend <= VOCEndDate[VOCList[j]])
                    {
                        foreach (string val in haploMaxPropDic.Keys)
                        {
                            if (haploMaxPropDic[val] >= 0.05)
                            {
                                string output = VOCList[j] + "\t" + Weeklist[i].Weekdate + "\t" + val;
                                output += "\t" + HaploDic[val].representingLineage;

                                if (Weeklist[i].WeekVOCHaploCountDic[VOCList[j]].ContainsKey(val))
                                    output += "\t" + Convert.ToString(Convert.ToDouble(Weeklist[i].WeekVOCHaploCountDic[VOCList[j]][val]) / Weeklist[i].WeekVOCCountDic[VOCList[j]]);
                                else
                                    output += "\t0";
                                writeProp.WriteLine(output);
                            }
                        }
                    }
                }
            
                
            }
            //write.Close();
            //writeProp.Close();
            //writeCase.Close();

        }
        static void CalculateMatrix(string fold, string country)//计算突变频率矩阵
        {
            int i, j, k;
            foreach(string val in HaploDic.Keys)
            {
                foreach (string mut in  HaploDic[val].SpikeMutDic.Keys)
                {
                    if (!SpikeMutList.Contains(mut)) SpikeMutList.Add(mut);
                }
            }

            for (j = 0; j < VOCList.Count(); j++)
            {
                StreamWriter write = new StreamWriter(fold + "/Matrix/VOCMatrix." + country + "." + VOCList[j].Replace('/','_') + ".tsv");
                string output = VOCList[j] + "\tDate\tIfWindow\tWeekCases\tAccumualteCases\tCasesProp";
                for (i = 0; i < SpikeMutList.Count(); i++)
                    output += "\t" + SpikeMutList[i];
                write.WriteLine(output);

                double accumulateCases = 0;
                for(i=0;i<Weeklist.Count();i++)
                {
                    output = VOCList[j] + "\t" + Weeklist[i].Weekdate;
                    if (Weeklist[i].Weekstart >= VOCStartDate[VOCList[j]] && Weeklist[i].Weekend <= VOCEndDate[VOCList[j]])
                        output += "\tY";
                    else
                        output += "\tN";
                    double weeklycases = Convert.ToDouble(Weeklist[i].WeekVOCCountDic[VOCList[j]]) / Weeklist[i].WeekSeqList.Count() * Weeklist[i].CasesSmoothed;
                    output += "\t" + Convert.ToString(weeklycases);

                    accumulateCases += weeklycases;
                    output += "\t" + Convert.ToString(accumulateCases);
                    output += "\t" + Convert.ToString(accumulateCases / VOCTotalCases[VOCList[j]]);

                    for (k = 0; k < SpikeMutList.Count(); k++)
                        if (Weeklist[i].WeekVOCCountDic[VOCList[j]] == 0)
                            output += "\t0";
                        else
                        {
                            int n=0;
                            foreach(string val in Weeklist[i].WeekVOCHaploCountDic[VOCList[j]].Keys)
                            {
                                if (HaploDic[val].SpikeMutDic.ContainsKey(SpikeMutList[k]))
                                    n += Weeklist[i].WeekVOCHaploCountDic[VOCList[j]][val];
                            }
                            output += "\t" + Convert.ToString(Convert.ToDouble(n) / Weeklist[i].WeekVOCCountDic[VOCList[j]]);
                        }
                    write.WriteLine(output);
                }
                write.Close();
            }
        }
        static void FindChangedMut(string fold, string country)//寻找并输出频率变化较大的突变
        {
            string[] filePaths = Directory.GetFiles("//NAS8500/g/VariationMutation/AiDMS/Spike0417/WeeklyMatrix/Matrix", "*", SearchOption.TopDirectoryOnly);
            int i, j, k;
            for(i=0;i<filePaths.Count();i++)
            {
                StreamReader read = new StreamReader(filePaths[i]);
                string line = read.ReadLine();
                string[] title = line.Split('\t');
                List<string> Mutlist = new List<string>();
                List<double> Maxfreq = new List<double>();
                List<double> Minfreq = new List<double>();
                List<double> T1freq = new List<double>();
                for(j=6;j<title.Length;j++)
                {
                    Mutlist.Add(title[j]);Maxfreq.Add(0);Minfreq.Add(1); T1freq.Add(0);
                }

                k = 0;
                line = read.ReadLine();
                while(line!=null)
                {
                    string[] line1 = line.Split('\t');
                    if (line1[2]=="Y")
                    {
                        k++;
                        for (j = 6; j < title.Length; j++)
                        {
                            if (Convert.ToDouble(line1[j]) > Maxfreq[j - 6]) Maxfreq[j - 6] = Convert.ToDouble(line1[j]);
                            if (Convert.ToDouble(line1[j]) < Minfreq[j - 6]) Minfreq[j - 6] = Convert.ToDouble(line1[j]);
                            if(k==1)T1freq[j-6]= Convert.ToDouble(line1[j]);
                        }
                    }
                    line = read.ReadLine();
                }

                //sort
                string tmps; double tmpd;
                for(j=0;j<Mutlist.Count();j++)
                    for(k=j+1;k<Mutlist.Count();k++)
                        if ((Maxfreq[j] - Minfreq[j]) < (Maxfreq[k] - Minfreq[k]))
                        {
                            tmps = Mutlist[j]; Mutlist[j] = Mutlist[k]; Mutlist[k] = tmps;
                            tmpd = Maxfreq[j]; Maxfreq[j] = Maxfreq[k]; Maxfreq[k] = tmpd;
                            tmpd = Minfreq[j]; Minfreq[j] = Minfreq[k]; Minfreq[k] = tmpd;
                            tmpd = T1freq[j]; T1freq[j] = T1freq[k]; T1freq[k] = tmpd;
                        }

                StreamWriter write = new StreamWriter(fold + "/FreqHeatmap." + title[0].Replace("/","_") + ".txt");
                StreamWriter writeanno = new StreamWriter(fold + "/FreqHeatmapChangeT1." + title[0].Replace("/", "_") + ".tsv");
                Dictionary<string, double> topMut = new Dictionary<string, double>();
                string outputHead = title[0] + "\t" + title[1] + "\t" + title[2] + "\t" + title[3] + "\t" + title[4] + "\t" + title[5];
                string output = "";
                k = 10;
                for (j = 0; j < k; j++)
                {
                    if (Mutlist[j]=="")
                    {
                        k++;
                        continue;
                    }
                    topMut.Add(Mutlist[j], T1freq[j]);
                }
                for (j = 6; j < title.Count(); j++)
                {
                    if (topMut.ContainsKey(title[j]))
                    {
                        output += "\t" + title[j];
                    }
                }
                writeanno.WriteLine(title[1] + output);
                write.WriteLine(title[1] + output);

                read.Close();
                read = new StreamReader(filePaths[i]);
                line = read.ReadLine();
                line = read.ReadLine();
                while(line!=null)
                {
                    string[] line1 = line.Split('\t');
                    output = line1[1];
                    for (j = 6; j < line1.Count(); j++)
                        if (topMut.ContainsKey(title[j]))
                            output += "\t" + Convert.ToString(Convert.ToDouble(line1[j])- topMut[title[j]]);
                    if (line1[2] == "Y")
                        writeanno.WriteLine(output);

                    output = line1[1];
                    for (j = 6; j < line1.Count(); j++)
                        if (topMut.ContainsKey(title[j]))
                            output += "\t" + line1[j];
                    if (line1[2] == "Y")
                        write.WriteLine(output);

                    line = read.ReadLine();
                }
                write.WriteLine("TMP\t-1\t1\t0\t1\t0\t1\t0\t1\t0\t1");
                writeanno.WriteLine("TMP\t-1\t1\t0\t1\t0\t1\t0\t1\t0\t1");

                read.Close();
                write.Close();
                writeanno.Close();
            }
            return;
        }
        static void CalculateGenomeNUCmut(string fold, string country)//计算genome nuc 突变数随时间的变化
        {
            int i, j, k;
            StreamWriter write = new StreamWriter(fold + "/VOCWeeklyNUC." + country + ".tsv");
            write.WriteLine("VOC\tDate\tCases\tAverageNUCMutCount");
            for (i = 0; i < VOCList.Count(); i++)
            {
                for (j = 0; j < Weeklist.Count; j++)
                {
                    if (Weeklist[j].Weekstart >= VOCStartDate[VOCList[i]] && Weeklist[j].Weekend <= VOCEndDate[VOCList[i]])
                    {
                        string output = VOCList[i] + "\t" + Weeklist[j].Weekdate;
                        output += "\t" + Convert.ToString(Convert.ToDouble(Weeklist[j].WeekVOCCountDic[VOCList[i]]) / Weeklist[j].WeekSeqList.Count() * Weeklist[j].CasesSmoothed);
                        
                        output += "\t" + Convert.ToString(Weeklist[j].WeekVOCGenomeNUC[VOCList[i]] / Weeklist[j].WeekVOCGenomeCount[VOCList[i]]);
                        write.WriteLine(output);
                    }
                }

            }
            write.Close();
            return;
        }
        static void Main(string[] args)
        {
            string country = "World";
            ReadingLineagenote("//NAS8500/g/VariationMutation/AiDMS/Spike0417/Data/lineage_note20240418.txt");//Read the lineage annotation
            ReadinCases("//NAS8500/g/VariationMutation/AiDMS/Spike0417/Data/WHO-COVID-19-global-data.csv", country);//Read weekly cases number, build time window
            ReadinSeq("//NAS8500/g/VariationMutation/AiDMS/Spike0417/Data", country);//Read sequence，assinged to each time window
            CalculateRepresentingLineage();//Find the representive lineage for each haplotype
            CalculateVOCCases("//NAS8500/g/VariationMutation/AiDMS/Spike0417/VOCCases", country);//Calculate the weekly cases of each clade

            ReadinVOCTimeWindow("//NAS8500/g/VariationMutation/AiDMS/Spike0417/Data/VOCTimeWindow.txt");//Read the analysis window for each clade
            CalculateHaploDiversity("//NAS8500/g/VariationMutation/AiDMS/Spike0417/Diversity", country);//Calculate diversity over time
            CalculateVOCHaplotype("//NAS8500/g/VariationMutation/AiDMS/Spike0417/Diversity", country);//Calculate haplotype distribution for each clade
            CalculateMatrix("//NAS8500/g/VariationMutation/AiDMS/Spike0417/WeeklyMatrix", country);//Calculate mutation frequency heatmap
            FindChangedMut("//NAS8500/g/VariationMutation/AiDMS/Spike0417/WeeklyMatrix/Heatmap", country);//Find muations with big frequency change
            CalculateGenomeNUCmut("//NAS8500/g/VariationMutation/AiDMS/Spike0417/Diversity", country);//Calculate genome nuc mutation number over time
            return;
        }
    }
}
