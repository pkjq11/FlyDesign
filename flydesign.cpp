#include "flydesign.h"
#define M_PI 3.14159265358979323846

FlyDesign::FlyDesign()
{
}

double FlyDesign::D2R(double deg){
return deg / 180.0 * 3.14159265;
}

bool FlyDesign::Isfit()
{
    if((H>H_ky_high)||H<H_ky_low)
        return false;
    else return true;
}

bool FlyDesign::DesignLine(const int& csfs_,const double& fi_,const double& H_ky_low_,const double& H_ky_high_,const double& H_,
                           vector<vector<Point3d<double> > > &Air_line,vector<vector<Point3d<double> > > &Direct_line)
{
    //参数赋值
    csfs=csfs_;
    fi=fi_;
    H_ky_low=H_ky_low_;
    H_ky_high=H_ky_high_;
    H=H_;

    if(In_points.size()<1||!Isfit())
    {
        flag=false;
        return false;
    }

    L=H/tan(D2R(fi))*csfs;
    if(In_points.size()==1)
    {
        double theta=0;
        vector<Point3d<double> > air_temp;
        while(theta<=2*M_PI)
        {
            Point3d<double> P_air_line;
            P_air_line.y=In_points[0].y+(H/tan(D2R(fi)))*sin(theta)/111319.55;
            P_air_line.x=In_points[0].x+(H/tan(D2R(fi)))*cos(theta)/(111319.55*cos(D2R(P_air_line.y)));
            P_air_line.z=In_points[0].z+H;
            air_temp.push_back(P_air_line);
            theta=theta+M_PI/50;
        }
        Air_line.push_back(air_temp);
        Direct_line.push_back(air_temp);
    }
    else{
        //此次更新只考虑每条航线的首尾两点进行作图
        for(size_t i=1;i<In_points.size();i++)
        {
            vector<Point3d<double> > air_temp;
            vector<Point3d<double> > direct_temp;
            direct_temp.push_back(In_points[i-1]);
            direct_temp.push_back(In_points[i]);
            double deta_Y,deta_X;
            if(In_points[i].y==In_points[i-1].y&&In_points[i].x!=In_points[i-1].x)
                deta_Y=L/111319.55;
            else if(In_points[i].x==In_points[i-1].x&&In_points[i].y!=In_points[i-1].y)
                deta_Y=0;
            else if(In_points[i].x==In_points[i-1].x&&In_points[i].y==In_points[i-1].y)
            {
                flag=false;
                return false;
            }
            else
            {
                double tanAlpha = (In_points[i].y-In_points[i-1].y)/(In_points[i].x-In_points[i-1].x);
                deta_Y=abs((L*cos(atan(tanAlpha)))/111319.55)*csfs;
            }
            int xd=In_points[i-1].x<In_points[i].x?1:-1;
            int yd=In_points[i-1].y<In_points[i].y?1:-1;
            for(vector<Point3d<double> >::iterator iter=direct_temp.begin();iter!=direct_temp.end();iter++)
            {
                Point3d<double> P_air_line;
                P_air_line.y=(*iter).y+deta_Y*xd;
                if(In_points[i].y==In_points[i-1].y&&In_points[i].x!=In_points[i-1].x)
                    deta_X=0;
                else if(In_points[i].x==In_points[i-1].x&&In_points[i].y!=In_points[i-1].y)
                    deta_X=L/(111319.55*cos(D2R(P_air_line.y)));
                else
                {
                    double tanAlpha = (In_points[i].y-In_points[i-1].y)/(In_points[i].x-In_points[i-1].x);
                    deta_X=abs((L*sin(atan(tanAlpha)))/(111319.55*cos(D2R(P_air_line.y))))*csfs;
                }
                P_air_line.x=(*iter).x-deta_X*yd;
                P_air_line.z=H;
                air_temp.push_back(P_air_line);
            }
            Air_line.push_back(air_temp);
            Direct_line.push_back(direct_temp);
        }
    }
    flag=true;
    return true;
}
bool FlyDesign::OutCoverBeam(const double& sita_,vector<vector<Point3d<double> > > &Direct_line,vector<vector<Point3d<double> > > &Cover_beam_near,vector<vector<Point3d<double> > > &Cover_beam_far)
{
    if(!flag||In_points.size()<=1)
    {
        flag=false;
        return false;
    }
    sita=sita_;
    for(size_t i=1;i<In_points.size();i++)
    {
        double L1 = (H*(tan(D2R(90-fi))-tan(D2R(90-fi-0.5*sita))));
        double L2 = (H*(tan(D2R(90-fi+0.5*sita))-tan(D2R(90-fi))));
        double deta_Y1,deta_Y2,deta_X1,deta_X2;
        if(In_points[i].y==In_points[i-1].y&&In_points[i].x!=In_points[i-1].x)
        {
            deta_Y1=L1*csfs/111319.55;
            deta_Y2=L2*csfs/111319.55;
        }
        else if(In_points[i].x==In_points[i-1].x&&In_points[i].y!=In_points[i-1].y)
            deta_Y1=deta_Y2=0;
        else if(In_points[i].x==In_points[i-1].x&&In_points[i].y==In_points[i-1].y)
        {
            flag=false;
            return false;
        }
        else
        {
            double tanAlpha = (In_points[i].y-In_points[i-1].y)/(In_points[i].x-In_points[i-1].x);
            deta_Y1=abs((L1*cos(atan(tanAlpha)))/111319.55)*csfs;
            deta_Y2=abs((L2*cos(atan(tanAlpha)))/111319.55)*csfs;
        }
        int xd=In_points[i-1].x<In_points[i].x?1:-1;
        int yd=In_points[i-1].y<In_points[i].y?1:-1;
        vector<Point3d<double> > near_temp;
        vector<Point3d<double> > far_temp;
        for(size_t j=0;j<Direct_line[i-1].size();j++)
        {
            Point3d<double> P_cover_near;
            Point3d<double> P_cover_far;
            P_cover_near.y=Direct_line[i-1][j].y+xd*deta_Y1;
            P_cover_far.y=Direct_line[i-1][j].y-xd*deta_Y2;
            if(In_points[i].y==In_points[i-1].y&&In_points[i].x!=In_points[i-1].x)
                deta_X1=deta_X2=0;
            else if(In_points[i].x==In_points[i-1].x&&In_points[i].y!=In_points[i-1].y)
            {
                deta_X1=L1*csfs/(111319.55*cos(D2R(P_cover_near.y)));
                deta_X2=L2*csfs/(111319.55*cos(D2R(P_cover_far.y)));
            }
            else if(In_points[i].x==In_points[i-1].x&&In_points[i].y==In_points[i-1].y)
            {
                flag=false;
                return false;
            }
            else
            {
                double tanAlpha = (In_points[i].y-In_points[i-1].y)/(In_points[i].x-In_points[i-1].x);
                deta_X1=abs((L1*sin(atan(tanAlpha)))/(111319.55*cos(D2R(P_cover_near.y))))*csfs;
                deta_X2=abs((L2*sin(atan(tanAlpha)))/(111319.55*cos(D2R(P_cover_far.y))))*csfs;
            }
            P_cover_near.x=Direct_line[i-1][j].x-yd*deta_X1;
            P_cover_far.x=Direct_line[i-1][j].x+yd*deta_X2;
            P_cover_near.z=P_cover_far.z=Direct_line[i-1][j].z;
            near_temp.push_back(P_cover_near);
            far_temp.push_back(P_cover_far);
        }
        Cover_beam_near.push_back(near_temp);
        Cover_beam_far.push_back(far_temp);
    }
    flag=true;
    return true;
}

bool FlyDesign::OutSwath(const double& AD_start_,const double& AD_length_,const double& fs_,const double& tao_,
                         vector<vector<Point3d<double> > > &Direct_line,vector<vector<Point3d<double> > > &Swath_start,vector<vector<Point3d<double> > > &Swath_end)
{
    if(!flag||In_points.size()<=1)
    {
        flag=false;
        return false;
    }
    AD_start=AD_start_;
    AD_length=AD_length_;
    fs=fs_;
    tao=tao_;
    double W=((AD_length*1000/fs-tao)*300);
    deta_M=abs(H/tan(D2R(fi)))-sqrt(pow(AD_start*300,2)-H*H);
    if(deta_M<0||W<0||W<deta_M)
      {
        flag=false;
        return false;
      }
    double deta_N=deta_M-W;//波束覆盖另一侧的相反数
    for(size_t i=1;i<In_points.size();i++)
    {
        double deta_Y1,deta_Y2,deta_X1,deta_X2;
        if(In_points[i].y==In_points[i-1].y&&In_points[i].x!=In_points[i-1].x)
        {
            deta_Y1=deta_M*csfs/111319.55;
            deta_Y2=deta_N*csfs/111319.55;
        }
        else if(In_points[i].x==In_points[i-1].x&&In_points[i].y!=In_points[i-1].y)
            deta_Y1=deta_Y2=0;
        else if(In_points[i].x==In_points[i-1].x&&In_points[i].y==In_points[i-1].y)
        {
            flag=false;
            return false;
        }
        else
        {
            double tanAlpha = (In_points[i].y-In_points[i-1].y)/(In_points[i].x-In_points[i-1].x);
            deta_Y1=abs((deta_M*cos(atan(tanAlpha)))/111319.55)*csfs;
            deta_Y2=abs((deta_N*cos(atan(tanAlpha)))/111319.55)*csfs;
        }
        int xd=In_points[i-1].x<In_points[i].x?1:-1;
        int yd=In_points[i-1].y<In_points[i].y?1:-1;
        vector<Point3d<double> > start_temp;
        vector<Point3d<double> > end_temp;
        for(size_t j=0;j<Direct_line[i-1].size();j++)
        {
            Point3d<double> P_swath_start;
            Point3d<double> P_swath_end;
            P_swath_start.z=P_swath_end.z=Direct_line[i-1][j].z;
            P_swath_start.y=Direct_line[i-1][j].y+xd*deta_Y1;
            P_swath_end.y=P_swath_start.y-xd*deta_Y2;
            if(In_points[i].y==In_points[i-1].y&&In_points[i].x!=In_points[i-1].x)
                deta_X1=deta_X2=0;
            else if(In_points[i].x==In_points[i-1].x&&In_points[i].y!=In_points[i-1].y)
            {
                deta_X1=deta_M*csfs/(111319.55*cos(D2R(P_swath_start.y)));
                deta_X2=deta_N*csfs/(111319.55*cos(D2R(P_swath_end.y)));
            }
            else if(In_points[i].x==In_points[i-1].x&&In_points[i].y==In_points[i-1].y)
            {
                flag=false;
                return false;
            }
            else
            {
                double tanAlpha = (In_points[i].y-In_points[i-1].y)/(In_points[i].x-In_points[i-1].x);
                deta_X1=abs((deta_M*sin(atan(tanAlpha)))/(111319.55*cos(D2R(P_swath_start.y))))*csfs;
                deta_X2=abs((deta_N*sin(atan(tanAlpha)))/(111319.55*cos(D2R(P_swath_end.y))))*csfs;
            }
            P_swath_start.x=Direct_line[i-1][j].x-yd*deta_X1;
            P_swath_end.x=P_swath_start.x+yd*deta_X2;
            start_temp.push_back(P_swath_start);
            end_temp.push_back(P_swath_end);
        }
        Swath_start.push_back(start_temp);
        Swath_end.push_back(end_temp);
    }
    flag=true;
    return true;

}

bool FlyDesign::OutNESZ(const double& LightVs_,const double& Fre_,const double& gain_,const double& Fn_,const double& Ls_,
                        const double& Va_,const double& NESZ_YQ_,const double& Ka_,const double& Kr_,const double& k_,
                        const double& T0_,const double& Res_,const double& TpRatio_,const double& Pfzh_,
                        vector<P_NESZ> &R_NESZ,vector<P_NESZ> &R_NESZ_YQ)
{
    if(!flag)
    {
        flag=false;
        return false;
    }
    LightVs=LightVs_;
    Fre=Fre_;
    gain=gain_;
    Fn=Fn_;
    Ls=Ls_;
    Va=Va_;
    NESZ_YQ=NESZ_YQ_;
    Ka=Ka_;
    Kr=Kr_;
    k=k_;
    T0=T0_;
    Res=Res_;
    TpRatio=TpRatio_;
    Pfzh=Pfzh_;
    double lamda=LightVs/Fre;
    double GtGr=pow(10,(gain/10))*pow(10,(gain/10));
    double temp1=2*pow(4*M_PI,3)*k*T0*pow(10,Fn/10)*pow(10,Ls/10)*(Va/3.6);
    double temp2=GtGr*pow(lamda,3)*Ka*Kr*Res*Pfzh*TpRatio;
    double R_min=sqrt((L*csfs-deta_M)*(L*csfs-deta_M)+H*H);
    double R_max=sqrt((L*csfs+deta_M)*(L*csfs+deta_M)+H*H);
    double R_temp=R_min;
    while(R_temp<=R_max)
    {
       double NESZ=10*log10(temp1*pow(R_temp,3)*cos(asin(H/R_temp))/temp2);
       R_NESZ.push_back(make_pair(R_temp,NESZ));
       R_NESZ_YQ.push_back(make_pair(R_temp,NESZ_YQ));
       R_temp=R_temp+1;
    }

    return true;

}

