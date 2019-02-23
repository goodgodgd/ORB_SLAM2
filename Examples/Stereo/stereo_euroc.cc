/**
* This file is part of ORB-SLAM2.
*
* Copyright (C) 2014-2016 Raúl Mur-Artal <raulmur at unizar dot es> (University of Zaragoza)
* For more information see <https://github.com/raulmur/ORB_SLAM2>
*
* ORB-SLAM2 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* ORB-SLAM2 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with ORB-SLAM2. If not, see <http://www.gnu.org/licenses/>.
*/


#include <iostream>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <chrono>

#include <opencv2/core/core.hpp>
#include <opencv2/ccalib/omnidir.hpp>

#include "System.h"
#include "Output.h"


using namespace std;

void LoadImages(const string &strPathLeft, const string &strPathRight, const string &strPathTimes,
                vector<string> &vstrImageLeft, vector<string> &vstrImageRight, vector<double> &vTimeStamps);

int main(int argc, char **argv)
{
    if(argc != 8)
    {
        cerr << endl << "Usage: ./stereo_euroc path_to_vocabulary path_to_settings path_to_left_folder path_to_right_folder path_to_times_file loop_closing_on output_file" << endl;
        return 1;
    }

    // Retrieve paths to images
    vector<string> vstrImageLeft;
    vector<string> vstrImageRight;
    vector<double> vTimeStamp;
    LoadImages(string(argv[3]), string(argv[4]), string(argv[5]), vstrImageLeft, vstrImageRight, vTimeStamp);
    
    ORB_SLAM2::Output::instance().set((atoi(argv[6]) != 0), argv[7]);

    if(vstrImageLeft.empty() || vstrImageRight.empty())
    {
        cerr << "ERROR: No images in provided path." << endl;
        return 1;
    }

    if(vstrImageLeft.size()!=vstrImageRight.size())
    {
        cerr << "ERROR: Different number of left and right images." << endl;
        return 1;
    }

    // Read rectification parameters
    cv::FileStorage fsSettings(argv[2], cv::FileStorage::READ);
    if(!fsSettings.isOpened())
    {
        cerr << "ERROR: Wrong path to settings" << endl;
        return -1;
    }

    cv::Mat K_l, K_r, P_l, P_r, R_l, R_r, D_l, D_r;
    fsSettings["LEFT.K"] >> K_l;
    fsSettings["RIGHT.K"] >> K_r;

    fsSettings["LEFT.P"] >> P_l;
    fsSettings["RIGHT.P"] >> P_r;

    fsSettings["LEFT.R"] >> R_l;
    fsSettings["RIGHT.R"] >> R_r;

    fsSettings["LEFT.D"] >> D_l;
    fsSettings["RIGHT.D"] >> D_r;

    int rows_l = fsSettings["LEFT.height"];
    int cols_l = fsSettings["LEFT.width"];
    int rows_r = fsSettings["RIGHT.height"];
    int cols_r = fsSettings["RIGHT.width"];

    ORB_SLAM2::Frame::mRowOffset = fsSettings["RowOffset"];
    float scaleTo8bit = float(fsSettings["ScaleTo8bit"]);
    cout << "Row offset between left and right: " << ORB_SLAM2::Frame::mRowOffset << endl;

    if(K_l.empty() || K_r.empty() || P_l.empty() || P_r.empty() || R_l.empty() || R_r.empty() || D_l.empty() || D_r.empty() ||
            rows_l==0 || rows_r==0 || cols_l==0 || cols_r==0)
    {
        cerr << "ERROR: Calibration parameters to rectify stereo are missing!" << endl;
        return -1;
    }

    int omni = fsSettings["Omni"];
    cv::Mat M1l,M2l,M1r,M2r;
    if(omni==1)
    {
        std::cout << "undistort omnidir camera" << std::endl;
        double xi_l = fsSettings["LEFT.xi"];
        double xi_r = fsSettings["RIGHT.xi"];
        cv::omnidir::initUndistortRectifyMap(K_l,D_l,xi_l,R_l,P_l,cv::Size(cols_l,rows_l),CV_32F,M1l,M2l,cv::omnidir::RECTIFY_PERSPECTIVE);
        cv::omnidir::initUndistortRectifyMap(K_r,D_r,xi_r,R_r,P_r,cv::Size(cols_r,rows_r),CV_32F,M1r,M2r,cv::omnidir::RECTIFY_PERSPECTIVE);
    }
    else
    {
        cv::initUndistortRectifyMap(K_l,D_l,R_l,P_l.rowRange(0,3).colRange(0,3),cv::Size(cols_l,rows_l),CV_32F,M1l,M2l);
        cv::initUndistortRectifyMap(K_r,D_r,R_r,P_r.rowRange(0,3).colRange(0,3),cv::Size(cols_r,rows_r),CV_32F,M1r,M2r);
    }

    const int nImages = vstrImageLeft.size();

    // Create SLAM system. It initializes all system threads and gets ready to process frames.
    ORB_SLAM2::System SLAM(argv[1],argv[2],ORB_SLAM2::System::STEREO,true);

    // Vector for tracking time statistics
    vector<double> vTimesTrack;
    vTimesTrack.resize(nImages);

    cout << endl << "-------" << endl;
    cout << "Start processing sequence ..." << endl;
    cout << "Images in the sequence: " << nImages << endl << endl;

    // Main loop
    cv::Mat imLeft, imRight, imLeftRect, imRightRect;
    // char filename[50];
    for(int ni=0; ni<nImages; ni++)
    {
        // Read left and right images from file
        imLeft = cv::imread(vstrImageLeft[ni],CV_LOAD_IMAGE_UNCHANGED);
        imRight = cv::imread(vstrImageRight[ni],CV_LOAD_IMAGE_UNCHANGED);

        if(scaleTo8bit > 1.f && imLeft.type()==2)
        {
            // cout << "16bit images are scaled by 1./" << scaleTo8bit << endl;
            imLeft.convertTo(imLeft, CV_8U, 1.f/scaleTo8bit);
            imRight.convertTo(imRight, CV_8U, 1.f/scaleTo8bit);
        }

        if(imLeft.empty())
        {
            cerr << endl << "Failed to load image at: "
                 << string(vstrImageLeft[ni]) << endl;
            return 1;
        }

        if(imRight.empty())
        {
            cerr << endl << "Failed to load image at: "
                 << string(vstrImageRight[ni]) << endl;
            return 1;
        }

#ifdef COMPILEDWITHC11
        std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
#else
        std::chrono::monotonic_clock::time_point t1 = std::chrono::monotonic_clock::now();
#endif

        cv::remap(imLeft,imLeftRect,M1l,M2l,cv::INTER_LINEAR);
        cv::remap(imRight,imRightRect,M1r,M2r,cv::INTER_LINEAR);
        // sprintf(filename, "/data/output/dump/L%05d.png", ni);
        // cv::imwrite(filename, imLeftRect);
        // sprintf(filename, "/data/output/dump/R%05d.png", ni);
        // cv::imwrite(filename, imRightRect);

        double tframe = vTimeStamp[ni];

        // Pass the images to the SLAM system
        SLAM.TrackStereo(imLeftRect,imRightRect,tframe);

#ifdef COMPILEDWITHC11
        std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
#else
        std::chrono::monotonic_clock::time_point t2 = std::chrono::monotonic_clock::now();
#endif

        double ttrack= std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(t2 - t1).count();
        vTimesTrack[ni]=ttrack;

#ifndef VO_BENCH_ON
        ttrack= std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
        vTimesTrack[ni]=ttrack;

        // Wait to load the next frame
        double T=0;
        if(ni<nImages-1)
            T = vTimestamps[ni+1]-tframe;
        else if(ni>0)
            T = tframe-vTimestamps[ni-1];

        if(ttrack<T)
            usleep((T-ttrack)*1e6);
#endif
    }

    // Stop all threads
    SLAM.Shutdown();

    // Tracking time statistics
    sort(vTimesTrack.begin(),vTimesTrack.end());
    double totaltime = 0;
    for(int ni=0; ni<nImages; ni++)
    {
        totaltime+=vTimesTrack[ni];
    }
    cout << "-------" << endl << endl;
    cout << "total frames: " << nImages << endl;
    cout << "mean tracking time: " << totaltime/nImages << "ms" << endl;

    // Save camera trajectory
    SLAM.SaveTrajectoryTUM(ORB_SLAM2::Output::instance().outfile, vTimesTrack);

    return 0;
}

void LoadImages(const string &strPathLeft, const string &strPathRight, const string &strPathTimes,
                vector<string> &vstrImageLeft, vector<string> &vstrImageRight, vector<double> &vTimeStamps)
{
    ifstream fTimes;
    fTimes.open(strPathTimes.c_str());
    if(!fTimes.is_open())
    {
    	cout << "cannot open time file: " << strPathTimes << endl;
    	throw 1;
    }

    vTimeStamps.reserve(5000);
    vstrImageLeft.reserve(5000);
    vstrImageRight.reserve(5000);
    while(!fTimes.eof())
    {
        string s;
        getline(fTimes,s);
        if(!s.empty())
        {
            stringstream ss;
            ss << s;
            vstrImageLeft.push_back(strPathLeft + "/" + ss.str() + ".png");
            vstrImageRight.push_back(strPathRight + "/" + ss.str() + ".png");
            double t;
            ss >> t;
            vTimeStamps.push_back(t/1e9);

        }
    }
}
