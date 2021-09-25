/*==============================================================================

Copyright 2017 Dalian University of Technology .
All rights reserved

================================================================================
-- Please append file description informations here --
PdProduct
1. Set calculate parameters for SOLVER
2. TPdNode can be : TPdNodeBBPD / TPdNodeNOSPD / TPdNodeBeamPD
3. TSolver can be : TSolverBBPD / TSolverNOSPD / TSolverBeamPD
================================================================================
Date            Name                    Description of Change
2017/12/24		Zheng Guojun			Create
2019/07/31		Zheng Guojun			Unified the solver call process
$HISTORY$
================================================================================
*/

#ifndef DLUT_SAE_PERIDYNAMIC_PD_PRODUCT_HXX_20171224
#define DLUT_SAE_PERIDYNAMIC_PD_PRODUCT_HXX_20171224

#include "pd_serialize.hxx"
#include "pd_solver_beampd.hxx"

#include <time.h>

namespace DLUT
{
	namespace SAE
	{
		namespace PERIDYNAMIC
		{
			template<typename TSolver>
			class PdProduct
			{
			public:
				PdProduct() {}
				~PdProduct() {}

			public:
				void		SetFilePath(string ip_lsdyna_info)
				{
					m_str_lsdyna_file = ip_lsdyna_info;
					string::size_type index = m_str_lsdyna_file.find_last_of('.');
					if (m_str_lsdyna_file.npos != index)
					{
						m_str_output_ascii = m_str_lsdyna_file.substr(0, index) + ".ascii";
					}
				}
				void		SetRatioOfHorizonMeshsize(double ratio_h_m)
				{
					DLUT::SAE::PERIDYNAMIC::RATIO_OF_HORIZON_MESHSIZE = ratio_h_m;
				}
				void		SetConstantHorizon(double horizon)
				{
					DLUT::SAE::PERIDYNAMIC::HORIZON = horizon;
				}
				void		SetUsedConstantHorizon(bool used)
				{
					DLUT::SAE::PERIDYNAMIC::USE_CONSTANT_HORIZON = used;
				}
				void		SetLoadStep(int load_step)
				{
					m_load_step = load_step;
					if (m_load_step < 1)
					{
						m_load_step = 1;
					}
				}
				void		SetMaxIteratorNums(int max_iter_nums)
				{
					DLUT::SAE::PERIDYNAMIC::MAX_INTERATOR_NUMS = max_iter_nums;
				}
				void		SetConvergenceFactor(double conv_factor)
				{
					DLUT::SAE::PERIDYNAMIC::CONVERGENCE_FACTOR = conv_factor;
				}
				void		SetTimeStep(double time_step)
				{
					m_time_step = time_step;
					if (m_time_step < 0)
					{
						m_time_step = 1;
					}
				}
				void		SetPlotFrames(int plot_frames)
				{
					m_plot_frames = plot_frames;
				}
				void		SetIteratorNums(int iterator_nums)
				{
					m_iterator_nums = iterator_nums;
				}
			public:
				void		ExplicitAnalysis()
				{
					int plot_frames = (int)(ceil((double)m_iterator_nums / m_plot_frames));
					//	Step 1: Read Pre-Process information
					fun_initPreData();
					seri_ascii.OutputHwAscii(0);
						
					//	Step 2: Refresh initial velocity/Acceleration/
					m_pd_model.RefreshInitVelocityInfo();									

					TSolver solver;
					solver.Attach(m_pd_model);
					//	Step 4: Solution
					double start = clock();	
					for (int tt = 1; tt <= m_iterator_nums; ++tt)
					{
						//	Step 3: Refresh body force
						m_pd_model.RefreshLoadNodeInfo(tt * m_time_step);

						//	Step 4.0: Load
						m_pd_model.RefreshPreMotionInfo(tt * m_time_step);

						//	Step 4.1: Calculate all the PD force 
						solver.ExplicitSolve(m_time_step);

						//	Step 4.3: Update boundary/motion information
						m_pd_model.RefreshBoundarySpcInfo();
						//	Step 4.4: Check fracture
						solver.RefreshFracture();

						//	Step 4.5: Output information (INTERVAL = tt*m_dt)
						if ((tt % plot_frames) == 0)
						{
							seri_ascii.OutputHwAscii(tt);
							cout << "Iterator numbers = " << setw(6) << tt;
							cout << ", and total calculate time = \t" << (clock() - start) / 1000. << endl;
						}
					}
					
					fun_finishOut();

//					cout << "Finished, and the total calculate time = \t" << (clock() - start) / 1000. << endl;
				}
				void		ImplicitAnalysis()
				{
					double start, end, cost;
					start = clock();

//					cout << "Static Analysis:" << endl;
					//	Step 1: Read Pre-Process information
					fun_initPreData();	
					seri_ascii.OutputHwAscii(0);
					TSolver solver;
					solver.Attach(m_pd_model);
					
					for (int cur_step = 1; cur_step <= m_load_step; ++cur_step)
					{
						double cur_step_start = clock();

						int ITERATOR_NUM = 0;
						bool is_convergenced = solver.ImplicitSolve(cur_step, m_load_step, ITERATOR_NUM);

						double cur_step_cost = (clock() - cur_step_start) / 1000;
						if (is_convergenced)
						{
							seri_ascii.OutputHwAscii(cur_step);

							printf("%s%5d%s%s%5d%s%5f\n", "Current step=", cur_step, " Convergenced, ", "ITERATOR NUMS=", ITERATOR_NUM, ", CALCULATE TIME=", cur_step_cost);
						}
						else
						{
							seri_ascii.OutputHwAscii(cur_step);

							printf("%s%5d%s%s%5d%s%5f\n", "Current step=", cur_step, " Not convergenced, ", "ITERATOR NUMS=", ITERATOR_NUM, ", CALCULATE TIME=", cur_step_cost);

							break;
						}
					}

					end = clock();
					cost = end - start;
					cout << endl;
					cout << "Finished, and the total time is " << cost / 1000. << endl;

					fun_finishOut();
				}
			private:
				void		fun_initPreData()
				{
//					cout << "Begin to read Pre-Process datas..." << endl;
					//	Read information from *KEY files (LSDYNA)
					m_pd_model.Initialize();
					SerializeLsdyna dyna_seri;
					dyna_seri.ReadLsdynaFile(m_str_lsdyna_file, m_pd_model);

					//	读入信息之后更新PART的信息
					m_pd_model.UpdatePartInfo();
					//	Update family information for all PD nodes
					m_pd_model.UpdateFamilyInParts();

					//	Update initial damage value for all PD nodes
					m_pd_model.GenerateInitDamage();
					//	Generate initial crevice and delete the corresponding family information
					m_pd_model.GenerateInitCrevice();
//					cout << "Finished to read Pre-Process datas..." << endl;

					//	Initialize output serialize.
					seri_ascii.AttachModel(m_pd_model);
					seri_ascii.AttachFile(m_str_output_ascii);

					string::size_type index = m_str_lsdyna_file.find_last_of('.');
					if (m_str_lsdyna_file.npos != index)
					{
						string lsdyna_file_new = m_str_lsdyna_file.substr(0, index) + "_discret.key";
						dyna_seri.WriteLsdynaFile(lsdyna_file_new, m_pd_model);
					}
				}
				void		fun_finishOut()
				{
					seri_ascii.DetachFile();
				}

			private:
				TPdModel	m_pd_model;				//	PD data of everything
			private:
				int			m_load_step;			//	Numbers of Load step
				double		m_time_step;			//	Time step
				int			m_iterator_nums;		//	Numbers of Iterator at every load step
				int			m_plot_frames;			//	Numbers of time interval for output

			private:
				string		m_str_lsdyna_file;		//	Filename of other information
				string		m_str_output_ascii;		//	Filename of output information
			private:
				SerializeHwAscii		seri_ascii;
			};
		} // END OF PERIDYNAMIC
	} // END OF DLUT
} // END OF SAE

#endif