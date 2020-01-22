//============================================================================
// Name        : BioM.cpp
// Author      : Sentieon
// Version     :
// Copyright   : 
// Description : C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include <map>
#include <set>
#include <vector>
#include <assert.h>
#include <algorithm>
#include <string>
#include <iomanip>
#include "lpsolve/lp_lib.h"
using namespace std;

vector<char*> tokenize(char* line, const char sep) {
    vector<char*> result;
    char* t = line;
    while(char* t1 = strchr(t, sep)) {
        *t1 = 0;
        result.push_back(t);
        t = t1+1;
    }
    result.push_back(t);
    return result;
}

double calcDist(const vector<double>& S1, const vector<double>& S2, const double* Coef) {
    double dist = 0.0;
    for(size_t i=0; i<S1.size(); i++) {
        const double delta = S1[i]-S2[i];
        dist += delta*delta*Coef[i];
    }
    return dist;
}

void trimspace(vector<vector<string> >& v) {
    for(size_t i = 0; i<v.size(); i++) {
        for(size_t j=0; j<v[i].size(); j++) {
            if(v[i][j]==" ") v[i][j]="";
        }
    }
}

typedef struct {
    int Surv_Stat;
    int UnProcessed;
} PATIENT_CLI_OUT;

typedef struct {
    size_t min_max, geneindex, cnt;
    double diff;
} BestRecord;

vector<vector<string> > ReadFile(const char* fname, vector<string>& firstcolumn, const char delimit, vector<string>& firstline, const int RDelta, const int mode) {
    char *line = NULL;
    size_t len = 0;
    ssize_t read = -1;
    vector<vector<string> > results;
    FILE* fp = fopen(fname, "r");
    size_t linecnt=0;

    firstcolumn.clear();
    while ((read = getline(&line, &len, fp)) != -1) {
        const unsigned int oldlen = strlen(line);
        assert(line[oldlen-1]=='\n');
        if(line[oldlen-2]=='\r') line[oldlen-2] = 0;
        else line[oldlen-1] = 0;

        vector<char*> tokens = tokenize(line, delimit);
        if(!linecnt) {
            firstline.clear();
            for(size_t i=RDelta; i<tokens.size(); i++) firstline.push_back(tokens[i]);
        }
        if(linecnt || mode==1) {
            firstcolumn.push_back(tokens[0]);
            vector<string> a;
            for(size_t i=1; i<tokens.size(); i++) a.push_back(tokens[i]);
            results.push_back(a);
            assert(a.size()==firstline.size());
        }
        linecnt++;
    }
    if (line)
        free(line);
    fclose(fp);
    cerr<<"sample count: "<<firstline.size()<<" vs. gene count: "<<firstcolumn.size()<<endl;
    return results;
}

int checkIndex(const vector<size_t>& t, const vector<PATIENT_CLI_OUT>& patient_cli_o, size_t& cnt) {
    int r = -1;
    cnt = 0;
    for(size_t i=0; i<t.size(); i++) {
        const PATIENT_CLI_OUT& pco = patient_cli_o[t[i]];
        if(pco.UnProcessed) {
            if(r<0) {
                r=pco.Surv_Stat;
                cnt=1;
            }
            else if(r==pco.Surv_Stat) {
                cnt++;
            }
            else {
                cnt=0;
                return -2;
            }
        }
    }
    return r;
}

double FindMin(double &mind, size_t& min_cnt, size_t& min_index, const map<double,vector<size_t> >& t, const vector<PATIENT_CLI_OUT>& patient_cli_o, const double margin) {
    assert(margin>=0.0);
    min_index = 0;
    mind = 0.0;
    min_cnt = 0;
    double MinThres = 0.0;
    map<double,vector<size_t> >::const_iterator it = t.begin();
    for(; it != t.end(); it++) {
        size_t cnt = 0;
        const int r = checkIndex(it->second, patient_cli_o, cnt);
        if(r==-2) {
            mind = it->first;
            return mind;
        }
        if(r>=0) {
            min_index = r;
            MinThres = it->first;
            break;
        }
    }
    for(it++; it != t.end(); it++) {
        size_t cnt = 0;
        const int r = checkIndex(it->second, patient_cli_o, cnt);
        if(r==-2 || (r>=0 && r!=(int)min_index)) {
            mind = it->first;
            break;
        }
    }
    for(it = t.begin(); it != t.end(); it++) {
        if(it->first + margin >= mind) return MinThres;
        size_t cnt = 0;
        const int r = checkIndex(it->second, patient_cli_o, cnt);
        assert(r==-1 || r==(int)min_index);
        min_cnt += cnt;
    }
    return MinThres;
}

double FindMax(double &maxd, size_t& max_cnt, size_t& max_index, const map<double,vector<size_t> >& t, const vector<PATIENT_CLI_OUT>& patient_cli_o, const double margin) {
    assert(margin>=0.0);
    max_index = 0;
    maxd = 0.0;
    max_cnt = 0;
    double MaxThres = 0.0;
    map<double,vector<size_t> >::const_reverse_iterator it = t.rbegin();
    for(; it != t.rend(); it++) {
        size_t cnt = 0;
        const int r = checkIndex(it->second, patient_cli_o, cnt);
        if(r==-2) {
            maxd = it->first;
            return maxd;
        }
        if(r>=0) {
            max_index = r;
            MaxThres = it->first;
            break;
        }
    }
    for(it++; it != t.rend(); it++) {
        size_t cnt = 0;
        const int r = checkIndex(it->second, patient_cli_o, cnt);
        if(r==-2 || (r>=0 && r!=(int)max_index)) {
            maxd = it->first;
            break;
        }
    }
    for(it = t.rbegin(); it != t.rend(); it++) {
        if(it->first <= margin + maxd) return MaxThres;
        size_t cnt = 0;
        const int r = checkIndex(it->second, patient_cli_o, cnt);
        assert(r==-1 || r==(int)max_index);
        max_cnt += cnt;
    }
    return MaxThres;
}

size_t FlipIndex(const vector<size_t>& t, vector<PATIENT_CLI_OUT>& patient_cli_o, const size_t index) {
    size_t cnt = 0;
    for(size_t i=0; i<t.size(); i++) {
        PATIENT_CLI_OUT& pco = patient_cli_o[t[i]];
        if(pco.UnProcessed) {
            assert(pco.Surv_Stat == (int)index);
            pco.UnProcessed = 0;
            cnt++;
        }
    }
    return cnt;
}

size_t FlipMin(double &mind, double& thres, size_t& min_index, const map<double,vector<size_t> >& t, vector<PATIENT_CLI_OUT>& patient_cli_o, const double margin) {
    assert(margin>=0.0);
    mind = 0.0;
    min_index = 0;
    thres = 0.0;
    size_t min_cnt = 0;
    map<double,vector<size_t> >::const_iterator it = t.begin();
    for(; it != t.end(); it++) {
        size_t cnt = 0;
        const int r = checkIndex(it->second, patient_cli_o, cnt);
        if(r==-2) return 0;
        if(r>=0) {
            min_index = r;
            mind = it->first;
            break;
        }
    }
    for(it++; it != t.end(); it++) {
        size_t cnt = 0;
        const int r = checkIndex(it->second, patient_cli_o, cnt);
        if(r==-2 || (r>=0 && r!=(int)min_index)) {
            thres = it->first - margin;
            break;
        }
    }
    for(it = t.begin(); it != t.end(); it++) {
        if(it->first >= thres) return min_cnt;
        size_t cnt = 0;
        const int r = checkIndex(it->second, patient_cli_o, cnt);
        assert(r==-1 || r==(int)min_index);
        min_cnt += cnt;
        FlipIndex(it->second, patient_cli_o, min_index);
    }
    return min_cnt;
}

size_t FlipMax(double &maxd, double& thres, size_t& max_index, const map<double,vector<size_t> >& t, vector<PATIENT_CLI_OUT>& patient_cli_o, const double margin) {
    assert(margin>=0.0);
    maxd = 0.0;
    max_index = 0;
    thres = 0.0;
    size_t max_cnt = 0;
    map<double,vector<size_t> >::const_reverse_iterator it = t.rbegin();
    for(; it != t.rend(); it++) {
        size_t cnt = 0;
        const int r = checkIndex(it->second, patient_cli_o, cnt);
        if(r==-2) return 0;
        if(r>=0) {
            max_index = r;
            maxd = it->first;
            break;
        }
    }
    for(it++; it != t.rend(); it++) {
        size_t cnt = 0;
        const int r = checkIndex(it->second, patient_cli_o, cnt);
        if(r==-2 || (r>=0 && r!=(int)max_index)) {
            thres = it->first + margin;
            break;
        }
    }
    for(it = t.rbegin(); it != t.rend(); it++) {
        if(it->first <= thres) return max_cnt;
        size_t cnt = 0;
        const int r = checkIndex(it->second, patient_cli_o, cnt);
        assert(r==-1 || r==(int)max_index);
        max_cnt += cnt;
        FlipIndex(it->second, patient_cli_o, max_index);
    }
    return max_cnt;
}

size_t CountUnprocessed(const vector<PATIENT_CLI_OUT>& patient_cli_o) {
    size_t result = 0;
    for(size_t i = 0; i<patient_cli_o.size(); i++) {
        const PATIENT_CLI_OUT& pco = patient_cli_o[i];
        if(pco.UnProcessed) {
            result++;
        }
    }
    return result;
}

map<string, size_t> Vec2Map(const vector<string>& v) {
    map<string, size_t> result;
    for(size_t i = 0; i < v.size(); i++) result[v[i]] = i;
    return result;
}

typedef struct {
    string gene;
    size_t Surv_Stat, min_max;
    double o_thres;
} TreeModel;

vector<TreeModel> ConvertTreeModel(const vector<string>& v, const vector<vector<string> >& vd, size_t& lastclass, double& lastThres, string& lastGene) {
    vector<TreeModel> result(v.size());
    if(!v.size()) {
        lastclass = 1;
        lastThres = 0.0;
        lastGene = "";
        return result;
    }
    for(size_t i = 0; i < v.size(); i++) {
        TreeModel tm;
        tm.gene = v[i];
        tm.Surv_Stat = atoi(vd[i][0].c_str());
        tm.min_max = atoi(vd[i][1].c_str());
        tm.o_thres = atof(vd[i][2].c_str());
        result[i] = tm;
    }
    lastclass = 1 - result[v.size()-1].Surv_Stat;
    lastThres = result[v.size()-1].o_thres;
    lastGene = result[v.size()-1].gene;
    return result;
}

void LoadFeatureSet(set<string>& Features, const char* Model, const char* FeatureCnt) {
    vector<string> SigGene, Dummy;
    const vector<vector<string> > SigGeneCnt = ReadFile(Model, SigGene, '\t', Dummy, 1, 1);
    const int FC = atoi(FeatureCnt);
    const size_t TotalFeatureSelected = FC >=0 ? min(SigGeneCnt.size(), (size_t)FC) : SigGeneCnt.size();
    for(size_t i=0; i<TotalFeatureSelected; i++) {
        Features.insert(SigGene[i]);
    }
}

set<string> LoadFeatureSet(set<string>& Features, const char* KNNModel, const char* KNNFeatureCnt, const char* TreeModel, const char* TreeFeatureCnt) {
    LoadFeatureSet(Features, KNNModel, KNNFeatureCnt);
    LoadFeatureSet(Features, TreeModel, TreeFeatureCnt);
    return Features;
}

int main(int argc, char*argv[]) {
    if(argc<3) {
        std::cerr<<"Wrong Usage! "<<std::endl;
        exit(0);
    }
    srand(time(NULL));
    const int mode = atoi(argv[1]);
    if(mode == 0) { //Decision tree training
        vector<string> samples, Status;
        const vector<vector<string> > SurvLL = ReadFile(argv[2], samples, '\t', Status, 1, 0);
        assert(Status.size()==1);

        vector<string> samples1, Pheno;
        vector<vector<string> > PhenoLL = ReadFile(argv[3], samples1, '\t', Pheno, 1, 0);
        trimspace(PhenoLL);
        assert(Pheno.size()==4);
        assert(samples1==samples);

        vector<string> samples2, Gene;
        const vector<vector<string> > GeneLL = ReadFile(argv[4], samples2, '\t', Gene, 1, 0);
        assert(samples2==samples);

        const size_t FeatureThres = atoi(argv[5]);
        const double ThresMargin = atof(argv[6]);
        const size_t SampleCnt = samples.size();
        map<string, size_t> CategoryCnt;
        vector<map<string, string> > CateStr;
        for(size_t i = 0; i< Pheno.size(); i++) {
            map<string, string> x;
            for(size_t j=0; j<SampleCnt; j++) {
                const string CatStr = Pheno[i]+"_"+PhenoLL[j][i];
                x[PhenoLL[j][i]] = CatStr;
                CategoryCnt[CatStr]++;
            }
            CateStr.push_back(x);
        }
        vector<map<string, size_t> > CategoryFea;
        const size_t GeneCnt1 = Gene.size();
        size_t tGeneCnt = GeneCnt1;
        vector<string> GeneNames = Gene;
        for(size_t i = 0; i< CateStr.size(); i++) {
            const map<string, string>& Cx = CateStr[i];
            map<string, size_t> x;
            for(map<string, string>::const_iterator it = Cx.begin(); it != Cx.end(); it++) {
                GeneNames.push_back(it->second);
                x[it->first] = tGeneCnt++;
             }
            CategoryFea.push_back(x);
        }


        vector<PATIENT_CLI_OUT> patient_cli_o(SampleCnt);
        for(size_t i = 0; i<SampleCnt; i++) {
            PATIENT_CLI_OUT& p = patient_cli_o[i];
            p.UnProcessed = 1;
            p.Surv_Stat = atoi(SurvLL[i][0].c_str());
        }
        vector<map<double, vector<size_t> > > GeneOrder(tGeneCnt);
        size_t cnt[2]={0,0};
        vector<size_t> TSamples;
        for(size_t i = 0; i<SampleCnt; i++) {
            TSamples.push_back(i);
            cnt[patient_cli_o[i].Surv_Stat]++;
        }
        for(size_t j=0; j<GeneCnt1; j++) {
            map<double,vector<size_t> >& t=GeneOrder[j];
            for(size_t i = 0; i<TSamples.size(); i++) {
                t[atof(GeneLL[TSamples[i]][j].c_str())].push_back(TSamples[i]);
            }
        }
        for(size_t j = 0; j< Pheno.size(); j++) {
            const map<string, size_t>& CF = CategoryFea[j];
            for(size_t i = 0; i<TSamples.size(); i++) {
                map<string, size_t>::const_iterator it = CF.find(PhenoLL[TSamples[i]][j]);
                assert(it != CF.end());
                for(map<string, size_t>::const_iterator it1 = CF.begin(); it1 != CF.end(); it1++) {
                    map<string, size_t>::const_iterator it2 = CategoryCnt.find(GeneNames[it1->second]);
                    if(it2==CategoryCnt.end()) cerr<<it1->first<<endl;
                    GeneOrder[it1->second][it1 == it ? (it2->second>20? 10.0 : 0.1) : 0.0].push_back(TSamples[i]);
                }
            }
        }
        cerr<<"Status counts: ("<<cnt[0]<<", "<<cnt[1]<<")"<<endl;
        BestRecord LastBR;
        for(size_t iter = 0; iter<800 && cnt[0] && cnt[1]; iter++) {
            BestRecord BR;
            BR.cnt = 0;
            BR.geneindex = 0;
            BR.min_max = 0;
            BR.diff = 0.0;
            for(size_t j=FeatureThres; j<GeneOrder.size(); j++) {
                const map<double,vector<size_t> >& t=GeneOrder[j];
                double mind;
                size_t count_min, index_min;
                FindMin(mind, count_min, index_min, t, patient_cli_o, ThresMargin);
                if(count_min>BR.cnt && j<GeneCnt1) {
                    BR.cnt = count_min;
                    BR.geneindex = j;
                    BR.min_max = 0;
                }
                double maxd;
                size_t count_max, index_max;
                FindMax(maxd, count_max, index_max, t, patient_cli_o, ThresMargin);
                if(count_max>BR.cnt) {
                    BR.cnt = count_max;
                    BR.geneindex = j;
                    BR.min_max = 1;
                }
            }
            if(!BR.cnt) break;
            cerr<<"iter "<<iter<<": ("<<cnt[0]<<", "<<cnt[1]<<") Gene: "<<GeneNames[BR.geneindex]<<"; Training reduc (";
            size_t status_index;
            double threshold, o_thres;
            const size_t Fcnt = BR.min_max ? FlipMax(threshold, o_thres, status_index, GeneOrder[BR.geneindex], patient_cli_o, ThresMargin): FlipMin(threshold, o_thres, status_index, GeneOrder[BR.geneindex], patient_cli_o, ThresMargin);
            assert(Fcnt == BR.cnt);
            status_index ? (cerr<<"0, "<<BR.cnt<<") ") : (cerr<<BR.cnt<<", 0) ");
            BR.min_max ? cerr<<"upper; " : cerr<<"lower; ";
            cerr<<threshold<<", "<<o_thres<<endl;
            cout<<GeneNames[BR.geneindex]<<"\t"<<status_index<<"\t"<<BR.min_max<<"\t"<<std::setprecision(14)<<o_thres<<endl;
            cnt[status_index]-=BR.cnt;
            LastBR=BR;

            assert(CountUnprocessed(patient_cli_o)==cnt[0]+cnt[1]);
        }
    }
    else if(mode == 1) { //tree Model apply
        vector<string> Model_Param, Gene_Model;
        vector<vector<string> > ModelStr = ReadFile(argv[2], Gene_Model, '\t', Model_Param, 1, 1);
        size_t lastclass = 0;
        double lastThres = 0.0;
        string lastGene;
        vector<TreeModel> model = ConvertTreeModel(Gene_Model, ModelStr, lastclass, lastThres, lastGene);

        vector<string> Ph_SC1, Patient_Ph_SC1;
        vector<vector<string> > SC1_Ph = ReadFile(argv[3], Patient_Ph_SC1, '\t', Ph_SC1, 1, 0);
        trimspace(SC1_Ph);

        vector<string> RE_SC1, Patient_RE_SC1;
        vector<vector<string> > SC1_RE = ReadFile(argv[4], Patient_RE_SC1, '\t', RE_SC1, 1, 0);
        assert(Patient_RE_SC1==Patient_Ph_SC1);
        const map<string, size_t> GeneMap = Vec2Map(RE_SC1);

        const size_t SampleCnt = Patient_RE_SC1.size();

        cout<<"PATIENTID\tSURVIVAL_STATUS\tDEPTH"<<endl;
        for(size_t i=0; i<SampleCnt; i++) {
            cout<<Patient_RE_SC1[i]<<"\t";
            if(!model.size()) {
                cout<<"1\t0"<<endl;
            }
            else {
                bool lastflag = true;
                map<string, size_t> PhenoMap;
                for(size_t j=0; j<Ph_SC1.size(); j++) {
                    PhenoMap[Ph_SC1[j]+"_"+SC1_Ph[i][j]]=j;
                }
                for(size_t j=0; j<model.size(); j++) {
                    const TreeModel& M = model[j];
                    map<string, size_t>::const_iterator it = GeneMap.find(M.gene);
                    double gdata = 0.0;
                    if(it != GeneMap.end()) {
                        gdata = atof(SC1_RE[i][it->second].c_str());
                    }
                    else {
                        gdata = PhenoMap.find(M.gene) == PhenoMap.end() ? 0 : 10;
                    }
                    if(M.min_max) {
                        if(gdata>M.o_thres+1e-6) {
                            cout<<M.Surv_Stat<<"\t"<<j;
                            lastflag = false;
                            break;
                        }
                    }
                    else {
                        if(gdata+1e-6<M.o_thres) {
                            cout<<M.Surv_Stat<<"\t"<<j;
                            lastflag = false;
                            break;
                        }
                    }
                }
                if(lastflag) {
                    cout<<lastclass<<"\t"<<model.size();
                }
                cout<<endl;
            }
        }
    }
    else if(mode == 2) { //feature selection for nearest neighbor from training data
        vector<string> samples, Status;
        const vector<vector<string> > SurvLL = ReadFile(argv[2], samples, '\t', Status, 1, 0);
        assert(Status.size()==1);

        vector<string> samples1, Pheno;
        vector<vector<string> > PhenoLL = ReadFile(argv[3], samples1, '\t', Pheno, 1, 0);
        trimspace(PhenoLL);
        assert(Pheno.size()==4);
        assert(samples1==samples);

        vector<string> samples2, GeneNames;
        const vector<vector<string> > GeneLL = ReadFile(argv[4], samples2, '\t', GeneNames, 1, 0);
        assert(samples2==samples);

        const size_t tGeneCnt = GeneNames.size();
        const size_t SampleCnt = samples.size();

        vector<vector<double> > GeneD;
        vector<size_t> SampleSurviveFlag;
        size_t SStatus = 0;
        const double weight = atof(argv[5]);
        const int AllFlag = atof(argv[6]);
        const size_t Iter = atoi(argv[7]);
        for(size_t j=0; j<SampleCnt; j++) {
            const vector<string>& PV = PhenoLL[j];
            if(AllFlag==1 || PV[1] != "") {
                const size_t S = SurvLL[j][0] == "0" ? 0 : 1;
                SStatus += S;
                SampleSurviveFlag.push_back(S);
                GeneD.push_back(vector<double>(tGeneCnt));
                vector<double>& GD = GeneD[GeneD.size()-1];
                for(size_t i=0; i<tGeneCnt; i++) {
                    GD[i] = atof(GeneLL[j][i].c_str());
                }
            }
        }
        assert(GeneD.size()==SampleSurviveFlag.size());
        FILE* fp = fopen(argv[8], "w");
        fprintf(fp, "%lf\n", (double)SStatus/SampleSurviveFlag.size());
        fclose(fp);
        double* coef = new double[tGeneCnt];
        double* DEL = new double[tGeneCnt];
        for(size_t i=0; i<tGeneCnt; i++) coef[i] = 1.0/tGeneCnt;
        int* colno =  new int[tGeneCnt];
        REAL* row = new REAL[tGeneCnt];
        for(size_t iter=0; iter<Iter; iter++) {
            vector<vector<map<double, vector<size_t> > > > dist_matT(SampleSurviveFlag.size());
            for(size_t k1 = 0; k1 < SampleSurviveFlag.size(); k1++) {
                dist_matT[k1].resize(2);
            }
            for(size_t k1 = 0; k1 < SampleSurviveFlag.size(); k1++) {
                const vector<double>& S1 = GeneD[k1];
                for(size_t k2 = k1+1; k2 < SampleSurviveFlag.size(); k2++) {
                    const double dist = calcDist(S1, GeneD[k2], coef);
                    dist_matT[k1][SampleSurviveFlag[k2]][dist].push_back(k2);
                    dist_matT[k2][SampleSurviveFlag[k1]][dist].push_back(k1);
                }
            }
            size_t Inc0 = 0, Inc1 = 0;
            double MisDist = 0.0;;
            memset(DEL, 0, tGeneCnt*sizeof(double));
            for(size_t k1 = 0; k1 < SampleSurviveFlag.size(); k1++) {
                const size_t SClose = dist_matT[k1][0].begin()->first < dist_matT[k1][1].begin()->first ? 0 : 1;
                if(SampleSurviveFlag[k1]) {
                    Inc1 += 1-SClose;
                }
                else {
                    Inc0 += SClose;
                }
                const double factor = SampleSurviveFlag[k1] ? 1 : -weight;
                MisDist += factor * (dist_matT[k1][0].begin()->first - dist_matT[k1][1].begin()->first);
                const vector<double>& S1 = GeneD[k1];
                const vector<double>& S20 = GeneD[dist_matT[k1][0].begin()->second[0]];
                const vector<double>& S21 = GeneD[dist_matT[k1][1].begin()->second[0]];
                for(size_t i = 0; i<tGeneCnt; i++) {
                    const double delta0 = S20[i] - S1[i];
                    const double delta1 = S21[i] - S1[i];
                    DEL[i] += (delta0*delta0-delta1*delta1) * factor;
                }
            }

            /* We will build the model row by row
                   So we start with creating a model with 0 rows and 2 columns */
            lprec *lp = make_lp(0, tGeneCnt);
            assert(lp);

            {
                /* let us name our variables. Not required, but can be useful for debugging */
                for(size_t i=0; i<tGeneCnt; i++) {
                    char GeneName[GeneNames[i].length()+1];
                    memcpy(GeneName, GeneNames[i].c_str(), GeneNames[i].length()+1);
                    set_col_name(lp, i+1, GeneName);
                }
            }

            {
                set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */

                /* construct first row the sum of coefs <= 1 */
                /* all variables in lpsolve package are set to non-negative by default */

                size_t j=0;
                for(size_t i=0; i<tGeneCnt; i++) {
                    colno[j] = i + 1;
                    row[j++] = 1.0;
                }

                /* add the row to lpsolve */
                assert(add_constraintex(lp, j, row, colno, LE, 1.0));
                assert(j==tGeneCnt);
            }
            {
                /* set an upper bound for each coefficient */
                const double UppThres = .1;
                for(size_t i=0; i<tGeneCnt; i++) {
                    colno[0] = i + 1;
                }
                for(size_t i=0; i<tGeneCnt; i++) {
                    memset(row, 0, sizeof(*row)*tGeneCnt);
                    row[i] = 1.0;

                    /* add the row to lpsolve */
                    assert(add_constraintex(lp, tGeneCnt, row, colno, LE, UppThres));
                }
            }

            {
                set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */

                /* set the objective function coefs . DEL */
                size_t j = 0;

                for(size_t i=0; i<tGeneCnt; i++) {
                    colno[j] = i + 1; /* first column */
                    row[j++] = DEL[i];
                }

                /* set the objective in lpsolve */
                assert(set_obj_fnex(lp, j, row, colno));
            }

            {
                /* set the object direction to maximize */
                set_maxim(lp);

                /* just out of curioucity, now show the model in lp format on screen */
                /* this only works if this is a console application. If not, use write_lp and a filename */
                //write_LP(lp, stderr);
                /* write_lp(lp, "model.lp"); */

                /* I only want to see important messages on screen while solving */
                set_verbose(lp, IMPORTANT);

                /* Now let lpsolve calculate a solution */
                assert(solve(lp) == OPTIMAL);
            }

            {
                /* a solution is calculated, now lets get some results */
                /* variable values */
                get_variables(lp, row);
                const double stepsize = .1;
                for(size_t j = 0; j < tGeneCnt; j++) {
                    //printf("%s: %f\t", get_col_name(lp, j+1), row[j]);
                    coef[j] += (row[j] - coef[j]) * stepsize;
                }
                /* we are done now */
            }
            cerr<<"iter: "<<iter<<": ("<<Inc0<<","<<Inc1<<"); before "<<MisDist<<"; after "<<get_objective(lp);
            cerr<<endl;
            if(lp != NULL) {
                /* clean up such that all used memory by lpsolve is freed */
                delete_lp(lp);
            }
        }
        delete[] row;
        delete[] colno;
        for(size_t i=0; i<tGeneCnt; i++) {
            cout<<GeneNames[i]<<"\t"<<coef[i]<<endl;
        }
        delete[] DEL;
        delete[] coef;
    }
    else if(mode == 3) { // randomized KNN
        //Training data
        vector<string> samples, Status;
        const vector<vector<string> > SurvLL = ReadFile(argv[2], samples, '\t', Status, 1, 0);
        assert(Status.size()==1);

        vector<string> samples1, Pheno;
        vector<vector<string> > PhenoLL = ReadFile(argv[3], samples1, '\t', Pheno, 1, 0);
        trimspace(PhenoLL);
        assert(Pheno.size()==4);
        assert(samples1==samples);

        vector<string> samples2, Gene;
        const vector<vector<string> > GeneLL = ReadFile(argv[4], samples2, '\t', Gene, 1, 0);
        assert(samples2==samples);

        //Testing data
        vector<string> samples_t, Pheno_t;
        vector<vector<string> > PhenoL_t = ReadFile(argv[6], samples_t, '\t', Pheno_t, 1, 0);
        trimspace(PhenoL_t);
        assert(Pheno_t.size()==4);
        assert(Pheno==Pheno_t);

        vector<string> samples2_t, Gene_t;
        const vector<vector<string> > GeneL_t = ReadFile(argv[7], samples2_t, '\t', Gene_t, 1, 0);
        assert(samples2_t==samples_t);
        const map<string, size_t> GeneMap_t = Vec2Map(Gene_t);

        const size_t tGeneCnt = Gene.size();
        const size_t SampleCntR = samples.size();

        size_t TSS_CNT = 0, TTS = 0;
        vector<vector<double> > GeneD;
        const int AllFlag = atof(argv[8]);
        vector<size_t> SampleSurviveFlag;
        for(size_t j=0; j<SampleCntR; j++) {//Train
            const vector<string>& PV = PhenoLL[j];
            if(AllFlag==1 || PV[1] != "") {
                const size_t S = SurvLL[j][0] == "0" ? 0 : 1;
                SampleSurviveFlag.push_back(S);
                TTS+=S;
                GeneD.push_back(vector<double>(tGeneCnt));
                vector<double>& GD = GeneD[GeneD.size()-1];
                for(size_t i=0; i<tGeneCnt; i++) {
                    GD[i] = atof(GeneLL[j][i].c_str());
                }
            }
        }
        TSS_CNT = SampleSurviveFlag.size();
        cerr<<"Process "<<TTS<<" out of "<<SampleSurviveFlag.size()<<endl;
        assert(GeneD.size()==SampleSurviveFlag.size());
        vector<vector<double> > GeneD_Test;
        for(size_t j=0; j<samples2_t.size(); j++) {//Test
            const vector<string>& PV = PhenoL_t[j];
            if(PV[1] != "") {
                GeneD_Test.push_back(vector<double>(tGeneCnt));
                vector<double>& GD = GeneD_Test[GeneD_Test.size()-1];
                for(size_t i=0; i<tGeneCnt; i++) {
                    const map<string, size_t>::const_iterator Git = GeneMap_t.find(Gene[i]);
                    assert(Git != GeneMap_t.end());
                    GD[i] = atof(GeneL_t[j][Git->second].c_str());
                }
            }
        }
        cerr<<"Process "<<GeneD_Test.size()<<" testing samples"<<endl;
        vector<string> SigGene, Dummy;
        const vector<vector<string> > SigGeneCnt = ReadFile(argv[5], SigGene, '\t', Dummy, 1, 1);
        assert(Dummy.size()==1);

        const size_t TotalFeatureSelected = min(SigGeneCnt.size(), (size_t)atoi(argv[19]));
        int SignificantIndex[TotalFeatureSelected];
        double Coef[tGeneCnt];
        memset(Coef, 0, sizeof(double)*tGeneCnt);
        {
            const map<string, size_t> FeatureStrMap = Vec2Map(Gene);
            for(size_t i=0; i<TotalFeatureSelected; i++) {
                map<string, size_t>::const_iterator Fit = FeatureStrMap.find(SigGene[i]);
                assert(Fit != FeatureStrMap.end());
                SignificantIndex[i] = Fit->second;
                Coef[SignificantIndex[i]] = atof(SigGeneCnt[i][0].c_str());
            }
        }
        const size_t FeatureSize= atoi(argv[9]);
        assert(FeatureSize <= TotalFeatureSelected);

        const int NeighborCountMin = atoi(argv[10]);
        const double RatioMin = atof(argv[11]);

        const int NeighborCountMax = atoi(argv[12]);
        const double RatioMax = atof(argv[13]);

        const int NeighborCountDelta = atoi(argv[14]);
        const double RatioDelta = atof(argv[15]);

        const int DIST_SCALE_FLAG = atoi(argv[16]);
        const int RANDOMCOUNT = atoi(argv[17]);
        const size_t TOPFIXED = atoi(argv[18]);
        assert(TOPFIXED<=FeatureSize);
        srand(1);

        vector<pair<double, double> > MislabelVote(GeneD_Test.size());
        for(int RandC = 0; RandC < RANDOMCOUNT; RandC++) {
            vector<map<double, vector<size_t> > > dist_matT(TSS_CNT + GeneD_Test.size());
            map<int, int> map_array;
            for(size_t kd = 0; kd < TOPFIXED; kd++) map_array[SignificantIndex[kd]] = kd;
            const size_t SelectRange = TotalFeatureSelected-TOPFIXED;
            int tempindex[SelectRange];
            for(size_t kd = 0; kd < SelectRange; kd++) tempindex[kd] = TOPFIXED+kd;
            for(size_t kl = 0; kl < FeatureSize-TOPFIXED; kl++) {
                const int kd = int((double)rand()*(SelectRange-kl)/(RAND_MAX+1.0));
                assert(kd+kl<SelectRange);
                int temp = tempindex[kd];
                map_array[SignificantIndex[temp]] = temp;
                tempindex[kd] = tempindex[SelectRange-kl-1];
            }
            assert(map_array.size()==FeatureSize);
            double TCoef[tGeneCnt];
            memset(TCoef, 0, tGeneCnt*sizeof(double));
            for(map<int, int>::const_iterator it = map_array.begin(); it != map_array.end(); it++) {
                TCoef[it->first] = DIST_SCALE_FLAG ? Coef[it->first] : 1.0;
            }
            for(size_t k1 = 0; k1 < TSS_CNT; k1++) {
                const vector<double>& GK1 = GeneD[k1];
                for(size_t k2 = k1 + 1; k2 < TSS_CNT; k2++) {
                    const double dist_index = calcDist(GK1, GeneD[k2], TCoef);
                    dist_matT[k1][dist_index].push_back(k2);
                    dist_matT[k2][dist_index].push_back(k1);
                }
            }
            for(size_t k1 = 0; k1 < GeneD_Test.size(); k1++) {
                const vector<double>& GK1 = GeneD_Test[k1];
                for(size_t k2 = 0; k2 < TSS_CNT; k2++) {
                    const double dist_index = calcDist(GK1, GeneD[k2], TCoef);
                    dist_matT[k1+TSS_CNT][dist_index].push_back(k2);
                }
            }
            cerr<<RandC<<" out of "<<RANDOMCOUNT<<" with accuracy of: ";
            for(int NeighborCount=NeighborCountMin; NeighborCount<=NeighborCountMax; NeighborCount+=NeighborCountDelta) {
                for(double Ratio = RatioMin; Ratio <= RatioMax; Ratio += RatioDelta) {
                    vector<pair<double, double> > MislabelVoteScore(dist_matT.size());
                    for(size_t k = 0; k < dist_matT.size(); k++) {
                        vector<pair<size_t, double> > Neighbors;
                        double VoteCat[2] = {0.0, 0.0};
                        for(map<double, vector<size_t> >::const_iterator it1 = dist_matT[k].begin(); it1 != dist_matT[k].end() && (int)Neighbors.size()<NeighborCount; it1++) {
                            for(size_t l=0; l<it1->second.size() && (int)Neighbors.size()<NeighborCount; l++) {
                                VoteCat[SampleSurviveFlag[it1->second[l]]] += it1->first < .2 ? 5.0 : 1.0/it1->first;
                                Neighbors.push_back(make_pair(it1->second[l], it1->first));
                            }
                        }
                        MislabelVoteScore[k].first=Ratio * VoteCat[0];
                        MislabelVoteScore[k].second=VoteCat[1];
                    }
                    vector<pair<double, double> > MislabelVoteScoreF(dist_matT.size());
                    map<double, vector<size_t> > ScoreRank;

                    for(size_t k = 0; k < MislabelVoteScore.size(); k++) {
                        const pair<double, double>& MI = MislabelVoteScore[k];
                        const double TotalCL = MI.first + MI.second;
                        MislabelVoteScoreF[k].first = MI.first/TotalCL;
                        MislabelVoteScoreF[k].second = MI.second/TotalCL;
                        if(k < TSS_CNT)
                            ScoreRank[MI.second/TotalCL].push_back(k);
                    }
                    size_t TopMatched = 0, id = 0;
                    for(map<double, vector<size_t> >::const_reverse_iterator Sit = ScoreRank.rbegin(); Sit != ScoreRank.rend() && id < TTS; Sit++) {
                        for(size_t kkk=0; kkk<Sit->second.size() && id < TTS; kkk++, id++) {
                            if(SampleSurviveFlag[Sit->second[kkk]]) {
                                TopMatched++;
                            }
                        }
                    }
                    const double VotingPower = double(TopMatched)/TTS;
                    cerr<<VotingPower<<","<<NeighborCount<<","<<Ratio<<");";
                    for(size_t k = TSS_CNT; k < MislabelVoteScoreF.size(); k++) {
                        MislabelVote[k-TSS_CNT].first += MislabelVoteScoreF[k].first*VotingPower;
                        MislabelVote[k-TSS_CNT].second += MislabelVoteScoreF[k].second*VotingPower;
                    }
                }
            }
            cerr<<endl;
        }
        for(size_t j=0, k=0; j<samples2_t.size(); j++) {//Test
            const vector<string>& PV = PhenoL_t[j];
            cout<<samples2_t[j]<<"\t";
            if(AllFlag==1 || PV[1] != "") {
                cout<<MislabelVote[k].second/(MislabelVote[k].first+MislabelVote[k].second);
                k++;
            }
            else {
                cout<<"1";
            }
            cout<<endl;
        }
    }
    else if(mode == 4) { //ensemble of tree results and KNN results
        vector<string> samples, Dummy;
        vector<vector<string> > TreeRST = ReadFile(argv[2], samples, '\t', Dummy, 1, 0);
        assert(Dummy.size()==2);

        vector<string> samples1, Dummy1;
        vector<vector<string> > KnnRST = ReadFile(argv[3], samples1, '\t', Dummy1, 1, 1);
        assert(Dummy1.size()==1);
        assert(samples1 == samples);

        const size_t RankThres = atoi(argv[4]);
        vector<map<int, map<double, vector<size_t> > > > TreeRST_Sort(2);
        for(size_t i=0; i<TreeRST.size(); i++) {
            const vector<string>& O=TreeRST[i];
            assert(O[0]=="1"||O[0]=="0");
            const int P = atoi(O[0].c_str());
            const int R = min(atoi(O[1].c_str()), (int)RankThres);
            TreeRST_Sort[P][R][atof(KnnRST[i][0].c_str())].push_back(i);
        }
        size_t Rank = 0;
        vector<size_t> FinalRank(samples.size());
        for(map<int, map<double, vector<size_t> > >::const_iterator it = TreeRST_Sort[0].begin(); it != TreeRST_Sort[0].end(); it++) {
            for(map<double, vector<size_t> >::const_iterator it1 = it->second.begin(); it1 != it->second.end(); it1++) {
                for(size_t i = 0; i < it1->second.size(); i++) {
                    FinalRank[it1->second[i]] = Rank++;
                }
            }
        }
        for(map<int, map<double, vector<size_t> > >::const_reverse_iterator it = TreeRST_Sort[1].rbegin(); it != TreeRST_Sort[1].rend(); it++) {
            for(map<double, vector<size_t> >::const_iterator it1 = it->second.begin(); it1 != it->second.end(); it1++) {
                for(size_t i = 0; i < it1->second.size(); i++) {
                    FinalRank[it1->second[i]] = Rank++;
                }
            }
        }
        for(size_t i=0; i<samples1.size(); i++) {
            cout<<samples1[i]<<"\t"<<((double)FinalRank[i]/(samples.size()-1))<<endl;
        }
    }
    else if(mode == 5) { // Compute performance
        vector<string> samples, Dummy;
        const vector<vector<string> > SurvLL = ReadFile(argv[2], samples, '\t', Dummy, 1, 0);
        assert(Dummy.size()==1);
        size_t SurVS = 0;
        vector<size_t> SampleSurviveFlag(samples.size());
        for(size_t j=0; j<samples.size(); j++) {
             const size_t S = SurvLL[j][0] == "0" ? 0 : 1;
             SampleSurviveFlag[j] = S;
             SurVS += S;
        }
        cerr<<"Load "<<samples.size()<<"("<<SurVS<<")"<<endl;
        vector<string> samples1;
        const vector<vector<string> > SurvLL1 = ReadFile(argv[3], samples1, '\t', Dummy, 1, 1);
        assert(Dummy.size()==1);
        assert(samples1==samples);
        map<double, vector<size_t> > ScoreRank;
        for(size_t i=0; i<samples1.size(); i++) {
            ScoreRank[atof(SurvLL1[i][0].c_str())].push_back(i);
        }
        const size_t T0 = (1-atof(argv[4])) * samples1.size()+.5;
        size_t W0 = 0, TP=0, FP=0;
        double AUC = 0.0;
        size_t WT0 = 0;
        for(map<double, vector<size_t> >::const_iterator it = ScoreRank.begin(); it != ScoreRank.end(); it++) {
            size_t DTP = 0;
            for(size_t i = 0; i<it->second.size(); i++, WT0++) {
                if(WT0<T0) W0+=(1-SampleSurviveFlag[it->second[i]]);
                DTP += 1-SampleSurviveFlag[it->second[i]];
            }
            AUC += (TP+TP+DTP)*(it->second.size()-DTP);
            TP += DTP;
            FP += it->second.size()-DTP;
        }
        const size_t T1 = samples1.size() - T0;
        AUC *= .5/SurVS/(samples1.size()-SurVS);
        const int TP1 = SurVS - T0 + W0;
        const double P1 = (double)TP1/T1;
        const double S1 = (double)TP1/SurVS;
        const double F11 = 2.0*TP1/(T1+SurVS);

        const double TP0 = W0;
        const double P0 = (double)TP0/T0;
        const double S0 = (double)TP0/(samples.size()-SurVS);
        const double F10 = 2.0*TP0/(T0 + samples.size()-SurVS);

        const double Accuracy = (double)(TP0+TP1)/samples.size();
        cout<<P0<<"\t"<<S0<<"\t"<<F10<<"\t"<<P1<<"\t"<<S1<<"\t"<<F11<<"\t"<<Accuracy<<"\t"<<AUC<<"\t"<<((double)TP1*TP0/(T1-TP1)/(T0-TP0))<<"\t[("<<TP0<<"\t"<<(T1-TP1)<<")\t("<<(samples.size()-T1-TP0)<<"\t"<<TP1<<")]"<<endl;
    }
    else if(mode == 6) { // generate final prediction
        vector<string> samples1, Dummy;
        const vector<vector<string> > SurvLL1 = ReadFile(argv[2], samples1, '\t', Dummy, 1, 1);
        map<double, vector<size_t> > ScoreRank;
        for(size_t i=0; i<samples1.size(); i++) {
            ScoreRank[atof(SurvLL1[i][0].c_str())].push_back(i);
        }
        const size_t T0 = (1-atof(argv[3])) * samples1.size()+.5;
        size_t WT0 = 0;
        vector<int> Predict(samples1.size());
        for(map<double, vector<size_t> >::const_iterator it = ScoreRank.begin(); it != ScoreRank.end(); it++) {
            for(size_t i = 0; i<it->second.size(); i++, WT0++) {
                if(WT0>=T0) Predict[it->second[i]] = 1;
            }
        }
        cout<<"PATIENTID\tSURVIVAL_STATUS"<<endl;
        for(size_t i=0; i<samples1.size(); i++) cout<<samples1[i]<<"\t"<<Predict[i]<<endl;
    }
    else if(mode == 7) { // Compute average performance
        const size_t Folder = atoi(argv[3]);
        vector<double> perf;
        for(size_t i = 0; i < Folder; i++) {
            char message[256];
            vector<string> Dummy1, Dummy2;
            snprintf(message, 256, "%s_%zu_perf", argv[2], i);
            const vector<vector<string> > SurvLL = ReadFile(message, Dummy1, '\t', Dummy2, 1, 1);
            assert(Dummy1.size()==1);
            assert(Dummy2.size()==12);
            if(!i) perf.resize(Dummy2.size()+1);
            perf[0]+=atof(Dummy1[0].c_str());
            size_t j=0;
            for(; j<Dummy2.size()-4; j++) perf[j+1] += atof(Dummy2[j].c_str());
            string& T = Dummy2[j++];
            assert(T[0]=='['&&T[1]=='(');
            perf[j] += atof(T.c_str()+2);
            T = Dummy2[j++];
            size_t TLen = T.length();
            assert(T[TLen-1]==')');
            T[TLen-1]=0;
            perf[j] += atof(T.c_str());
            T = Dummy2[j++];
            TLen = T.length();
            assert(T[0]=='(');
            perf[j] += atof(T.c_str()+1);
            T = Dummy2[j++];
            TLen = T.length();
            assert(T[TLen-2]==')'&&T[TLen-1]==']');
            T[TLen-2]=0;
            perf[j] += atof(T.c_str());
        }
        size_t i=0;
        for(; i<perf.size()-4; i++) cout<<(perf[i]/Folder)<<"\t";
        cout<<"[("<<(perf[i]/Folder)<<"\t"<<(perf[i+1]/Folder)<<")\t("<<(perf[i+2]/Folder)<<"\t"<<(perf[i+3]/Folder)<<")]"<<endl;
    }
}
