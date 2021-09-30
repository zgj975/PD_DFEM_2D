#ifndef DLUT_SAE_PERIDYNAMIC_FEM_DATABASE_HXX_20181130
#define DLUT_SAE_PERIDYNAMIC_FEM_DATABASE_HXX_20181130

#include <vector>
#include <set>
using namespace std;

namespace DLUT
{
	namespace SAE
	{
		namespace PERIDYNAMIC
		{
			typedef Eigen::Vector3d TCoordinate;
			typedef Eigen::Matrix<double, 6, 1> TForce;
			typedef Eigen::Matrix<double, 6, 1> TDisplacement;
			typedef Eigen::Matrix<double, 6, 1> TVelocity;
			typedef Eigen::Matrix<double, 6, 1> TAcceleration;
			typedef Eigen::Matrix<double, 6, 1> TStress;
			typedef Eigen::Matrix<double, 6, 1> TStrain;
			typedef Eigen::Matrix<double, 4, 1> TStress_Bond;
			typedef Eigen::Matrix<double, 4, 1> TStrain_Bond;
			typedef Eigen::Matrix<double, 48, 48>	MATRTIX_SINGLE_STIFFNESS;
			typedef Eigen::Matrix<double, 48, 1>	MATRIX_PD_FORCE;
			const int IP_COUNT_2D = 4;	//	2D单元有四个积分点
			const int IP_COUNT_1D = 3;	//	1D单元有三个积分点
			const int DOF = 6;			//	节点有3个平动自由度+3个转动自由度
			//	2D 
			const double S_IP_2D[5] = {1.0 / sqrt(3.0), 1.0 / sqrt(3.0), -1.0 / sqrt(3.0), -1.0 / sqrt(3.0), 0};
			const double T_IP_2D[5] = {1.0 / sqrt(3.0), -1.0 / sqrt(3.0), 1.0 / sqrt(3.0), -1.0 / sqrt(3.0), 0};
			//	Bond
			const double S_IP_1D[3] = { (-sqrt(0.6) + 1) / 2.0, 0.5,(sqrt(0.6) + 1) / 2.0 };
			const double H_IP_1D[3] = { 5.0 / 18.0, 8.0 / 18.0, 5.0 / 18.0 };

			enum ANALYSIS_ELEMENT_TYPE {FEM_ELEMENT, PD_ELEMENT, MORPHING_ELEMENT};
			enum MESH_ELEMENT_TYPE {TRIANGLE_ELEMENT, QUADRANGLE_ELEMENT};

			//	*BOUNDARY_SPC_NODE
			class TBoundarySpcNode
			{
			public:
				TBoundarySpcNode(int dof = 0)
					: m_dof(dof)
				{ /*Do nothing*/				
				}
				~TBoundarySpcNode() {/*Do nothing*/ }
			public:
				int&		Dof() { return m_dof; }
				int			Dof() const { return m_dof; }
			private:
				int			m_dof;
			};
			//	*LOAD_NODE_POINT
			class TLoadNodePoint
			{
			public:
				TLoadNodePoint(int dof = -1, double sf = 0.0, int lcid = 0)
					: m_dof(dof), m_sf(sf), m_lcid(lcid)
				{ /*Do nothing*/
				}
				~TLoadNodePoint() { /*Do nothing*/ }
			public:
				int&		Dof() { return m_dof; }
				int			Dof() const { return m_dof; }

				int&		Lcid() { return m_lcid; }
				int			Lcid() const { return m_lcid; }

				double		Sf() const { return m_sf; }
				double&		Sf() { return m_sf; }
			private:
				int			m_dof;
				int			m_lcid;
				double		m_sf;
			};
			//	*INITIAL_VELOCITY_NODE
			class TInitialVelocityNode
			{
			public:
				TInitialVelocityNode()
				{
					m_velocity.setZero();
				}
				~TInitialVelocityNode() { /*Do nothing*/ }
			public:
				TVelocity& Velocity() { return m_velocity; }
				const TVelocity& Velocity() const { return m_velocity; }
			private:
				int			m_nid;
				TVelocity	m_velocity;
			};
			//	*BOUNDARY_PRESCRIBED_MOTION_NODE
			class TBoundaryPrescribedMotion
			{
			public:
				TBoundaryPrescribedMotion() { /*Do nothing*/ }
				~TBoundaryPrescribedMotion() { /*Do nothing*/ }
			public:
				int& Lcid() { return m_lcid; }
				int			Lcid() const { return m_lcid; }

				int& Dof() { return m_dof; }
				int			Dof() const { return m_dof; }

				double& Sf() { return m_sf; }
				double		Sf() const { return m_sf; }

				int& Vda() { return m_vda; }
				int			Vda() const { return m_vda; }

				string& MotionType() { return motion_type; }
				const string& MotionType() const { return motion_type; }
			private:
				int			m_lcid;
				int			m_dof;
				int			m_vda;
				double		m_sf;
				string		motion_type;
			};

			//	Node
			class TNodeBase
			{
			public:
				TNodeBase() : m_p_local_id(NULL)
				{
					m_coordinate.setZero();
					m_nid_global = -1;
					m_set_adj_element_ids.clear();

					m_displacement.setZero();
					m_iter_displacement.setZero();
					m_delta_displacement.setZero();
					m_velocity.setZero();
					m_acceleration.setZero();

					m_force_of_inner.setZero();
					m_force_of_outer.setZero();

					m_load.clear();
					m_init_velocity.clear();
					m_boundary_spc.clear();
					m_boundary_pre_motion.clear();
				}
				TNodeBase(const TCoordinate& pt, int golbal_nid = -1)
				{
					m_p_local_id = NULL;	
					
					m_coordinate = pt;
					m_nid_global = golbal_nid;
					m_set_adj_element_ids.clear();

					m_displacement.setZero();
					m_iter_displacement.setZero();
					m_delta_displacement.setZero();
					m_velocity.setZero();
					m_acceleration.setZero();

					m_force_of_inner.setZero();
					m_force_of_outer.setZero();

					m_load.clear();
					m_init_velocity.clear();
					m_boundary_spc.clear();
					m_boundary_pre_motion.clear();
				}
				void				Dispose()
				{
					if (m_p_local_id)
					{
						delete m_p_local_id;
						m_p_local_id = NULL;
					}
					m_set_adj_element_ids.clear();
				}
			public:
				TCoordinate&		Coordinate() { return m_coordinate; }
				const TCoordinate&	Coordinate() const { return m_coordinate; }

				int&				Id() 
				{
					if (m_p_local_id == NULL) m_p_local_id = new int(-1);
					return *m_p_local_id;
				}
				int					Id() const
				{
					if (m_p_local_id == NULL)
						return -1;
					else
						return *m_p_local_id;
				}
				int&				IdGlobal() { return m_nid_global; }
				int					IdGlobal() const { return m_nid_global; }

				void				InsertAdjElement(int *eid) { m_set_adj_element_ids.insert(eid); }
				void				DeleteAdjElement(int* eid) { m_set_adj_element_ids.erase(eid); }
				set<int>			AdjElementIds() const
				{
					set<int> res;
					for(int* eid : m_set_adj_element_ids)
					{
						res.insert(*eid);
					}
					return res;
				}
				int					AdjElementCount() const
				{
					return (int)(m_set_adj_element_ids.size());
				}

			public:
				TDisplacement&			Displacement() { return m_displacement; }
				const TDisplacement&	Displacement() const { return m_displacement; }

				TVelocity&				Velocity() { return m_velocity; }
				const TVelocity&		Velocity()	const { return m_velocity; }

				TAcceleration&			Acceleration() { return m_acceleration; }
				const TAcceleration&	Acceleration() const	{ return m_acceleration; }
		
				TDisplacement&			IncrementalDisplacement() { return m_delta_displacement; }
				const TDisplacement&	IncrementalDisplacement() const { return m_delta_displacement; }

				TDisplacement&			IteratorDisplacement() { return m_iter_displacement; }
				const TDisplacement&	IteratorDisplacement() const { return m_iter_displacement; }

			public:
				TForce&					InnerForce() { return m_force_of_inner; }
				const TForce&			InnerForce() const { return m_force_of_inner; }

				TForce&					OuterForce() { return m_force_of_outer; }
				const TForce&			OuterForce() const { return m_force_of_outer; }
			public:
				vector<TBoundarySpcNode>&					BoundarySpcNode() { return m_boundary_spc; }
				const vector<TBoundarySpcNode>&				BoundarySpcNode() const { return m_boundary_spc; }
				vector<TLoadNodePoint>&						LoadNodePoint() { return m_load; }
				const vector<TLoadNodePoint>&				LoadNodePoint() const { return m_load; }
				void										ScaleLoadNodePoints(double scale)
				{
					for (TLoadNodePoint& loadNodePt: m_load)
					{
						loadNodePt.Sf() *= scale;
					}
				}
				vector<TInitialVelocityNode>&				InitVelocity() { return m_init_velocity; }
				const vector<TInitialVelocityNode>&			InitVelocity() const { return m_init_velocity; }
				vector<TBoundaryPrescribedMotion>&			BoundaryPreMotion() { return m_boundary_pre_motion; }
				const vector<TBoundaryPrescribedMotion>&	BoundaryPreMotion() const { return m_boundary_pre_motion; }
			private:
				TCoordinate				m_coordinate;						//	Coordinate 
				set<int*>				m_set_adj_element_ids;				//	Adjoint element set
				int						*m_p_local_id;						//	Id of this node
				int						m_nid_global;						//	Global Id of this node

				TDisplacement			m_displacement;						//	Displacement
				TVelocity				m_velocity;							//	Velocity
				TAcceleration			m_acceleration;						//	Acceleration

				TDisplacement			m_delta_displacement;				//	Incremental Displacement of current incremental step
				TDisplacement			m_iter_displacement;				//	Iterator Displacement
			private:
				TForce					m_force_of_inner;					//	Inner force of this node, such as PD bond force, FEM node force...
				TForce					m_force_of_outer;					//	Outer force of this node, such as Body force, Interface force...
			private:
				vector<TBoundarySpcNode>			m_boundary_spc;
				vector<TLoadNodePoint>				m_load;
				vector<TInitialVelocityNode>		m_init_velocity;
				vector<TBoundaryPrescribedMotion>	m_boundary_pre_motion;
			};

			//	Integration Point
			template<typename TStrain, typename TStress>
			class TIntegrationPointTemp
			{
			public:
				TIntegrationPointTemp(int index = 0)
				{
					m_coordinate.setZero();
					m_iter_displacement.setZero();
					m_delta_displacement.setZero();

					m_stress_current.setZero();
					m_strain_current.setZero();

					m_stress_laststep.setZero();
					m_strain_laststep.setZero();

					m_index = index;
				}
			public:
				TCoordinate&				Coordinate() { return m_coordinate; }
				const TCoordinate&			Coordinate() const { return m_coordinate; }
				TDisplacement&				IteratorDisplacement() { return m_iter_displacement; }
				const TDisplacement&		IteratorDisplacement() const { return m_iter_displacement; }
				TDisplacement&				IncrementalDisplacement() { return m_delta_displacement; }
				const TDisplacement&		IncrementalDisplacement() const { return m_delta_displacement; }
			public:
				TStress&					StressCurrent() { return m_stress_current; }
				const TStress&				StressCurrent() const { return m_stress_current; }
				TStrain&					StrainCurrent() { return m_strain_current; }
				const TStrain&				StrainCurrent() const { return m_strain_current; }

				TStress&					StressLaststep() { return m_stress_laststep; }
				const TStress&				StressLaststep() const { return m_stress_laststep; }
				TStrain&					StrainLaststep() { return m_strain_laststep; }
				const TStrain&				StrainLaststep() const { return m_strain_laststep; }
			public:
				int							Index() const { return m_index; }
			private:
				TCoordinate					m_coordinate;
				TDisplacement				m_iter_displacement;	//	积分点处当前时刻当前迭代步的位移增量
				TDisplacement				m_delta_displacement;	//	积分点处当前时刻当前增量步的位移增量
			private:
				TStress						m_stress_current;		//	积分点处当前时刻的应力
				TStrain						m_strain_current;		//	积分点处当前时刻的应变
				TStress						m_stress_laststep;		//	积分点处上一时刻的应力
				TStrain						m_strain_laststep;		//	积分点处上一时刻的应变
			private:
				int							m_index;				//	积分点的序号
			};

			typedef TIntegrationPointTemp<TStrain, TStress> TIntegrationPoint;
			//	Element
			class TElementBase
			{
			public:
				TElementBase(vector<TNodeBase>& vecNode) 
					: m_p_local_id(NULL), ref_nodes(vecNode) 
				{
					m_elem_type = FEM_ELEMENT;
					for (int index = 0; index <= IP_COUNT_2D; ++index)
					{
						m_IP.push_back(TIntegrationPoint(index));
					}
				}
				void				Dispose()
				{
					if (m_p_local_id)
					{
						delete m_p_local_id;
						m_p_local_id = NULL;
					}
				}
				const TElementBase& operator=(const TElementBase& right)
				{
					m_vec_nids = right.m_vec_nids;
					m_part_id = right.m_part_id;
					m_p_local_id = right.m_p_local_id;
					*m_p_local_id = *(right.m_p_local_id);
					m_eid_global = right.m_eid_global;
					m_elem_type = right.m_elem_type;

					m_side_length = right.m_side_length;
					m_side_a = right.m_side_a;
					m_side_b = right.m_side_b;
					m_area = right.m_area;
					m_thickness = right.m_thickness;
					m_local_coord_system = right.m_local_coord_system;

					m_IP = right.m_IP;

					return *this;
				}
			public:
				void				InsertNode(int* nid) { m_vec_nids.push_back(nid); }
				vector<int>			NodeIds() const 
				{
					vector<int> res;
					for (int* pnid : m_vec_nids)
					{
						res.push_back(*pnid);
					}
					return res;
				}
				int					NodeId(int pos) const 
				{
					assert(pos >= 0 && pos < m_vec_nids.size());

					return *(m_vec_nids[pos]);
				}
				const TNodeBase&	Node(int pos) const
				{
					int nid = NodeId(pos);
					return ref_nodes[nid];
				}
				int					NodeCount() const
				{
					return (int)(m_vec_nids.size());
				}
				void				UpdateNodes(vector<int*> nids)
				{
					m_vec_nids.clear();
					m_vec_nids = nids;
					RefreshDatas();
				}

				int&				Id() 
				{
					if (m_p_local_id == NULL)
						m_p_local_id = new int(-1);

					return *m_p_local_id;
				}
				int					Id() const 
				{
					if (m_p_local_id == NULL)
						return -1;
					else
						return *m_p_local_id;
				}
				int&				IdGlobal() { return m_eid_global; }
				int					IdGlobal() const { return m_eid_global; }

				int&				PartId() { return m_part_id; }
				int					PartId() const { return m_part_id; }

				ANALYSIS_ELEMENT_TYPE&		AnalysisElementType() { return m_elem_type; }
				ANALYSIS_ELEMENT_TYPE		AnalysisElementType() const { return m_elem_type; }
				MESH_ELEMENT_TYPE	MeshElementType() const
				{
					MESH_ELEMENT_TYPE res;
					switch (m_vec_nids.size())
					{
					case 3:
						res = TRIANGLE_ELEMENT;
						break;
					case 4:
						res = QUADRANGLE_ELEMENT;
						break;
					default:
						break;
					}
					
					return res;
				}
			public:
				double				SideLength() const { return m_side_length; }
				double				SideA() const { return m_side_a; }
				double				SideB() const { return m_side_b; }
				double				Area() const { return m_area; }
				double&				Thickness() { return m_thickness; }
				double				Thickness() const { return m_thickness; }

			public:
				double				Jacobi(double s, double t) const
				{		
					const Matrix3d& T = m_local_coord_system;
					Eigen::Matrix<double, 4, 4> ST;
					vector<int> nids = NodeIds();
					
					Vector4d Xc;
					Vector4d Yc;
					for (int i = 0; i < 4; ++i)
					{
						const TCoordinate& cord_local = T * ref_nodes[nids[i]].Coordinate();
						Xc(i) = cord_local.x();
						Yc(i) = cord_local.y();
					}

					ST.setZero();
					ST(0, 1) = 1 - t;
					ST(0, 2) = t - s;
					ST(0, 3) = s - 1;

					ST(1, 0) = t - 1;
					ST(1, 2) = s + 1;
					ST(1, 3) = -s - t;

					ST(2, 0) = s - t;
					ST(2, 1) = -s - 1;
					ST(2, 3) = t + 1;

					ST(3, 0) = 1 - s;
					ST(3, 1) = s + t;
					ST(3, 2) = -t - 1;

					return abs(1.0 / 8.0 * Xc.transpose() * ST * Yc);
				}
				TCoordinate			CoordinateInElement(double s, double t) const
				{
					const Matrix3d& T = m_local_coord_system;
					vector<int> nids = NodeIds();
					
					Vector4d N;
					N(0) = (1 - s) * (1 - t) / 4.0;
					N(1) = (1 + s) * (1 - t) / 4.0;
					N(2) = (1 + s) * (1 + t) / 4.0;
					N(3) = (1 - s) * (1 + t) / 4.0;

					Vector4d X;
					Vector4d Y;
					Vector4d Z;
					for (int i = 0; i < 4; ++i)
					{
						const Vector3d& cord_local = T * ref_nodes[nids[i]].Coordinate().block(0, 0, 3, 1);
						X(i) = cord_local.x();
						Y(i) = cord_local.y();
						Z(i) = cord_local.z();
					}

					TCoordinate res_local;
					res_local.x() = N.transpose() * X;
					res_local.y() = N.transpose() * Y;
					res_local.z() = N.transpose() * Z;

					TCoordinate res_global = T.transpose() * res_local;

					return res_global;
				}
				TDisplacement		DisplaceInElement(double s, double t) const 
				{
					const Matrix3d& T = m_local_coord_system;
					vector<int> nids = NodeIds();
					if (MeshElementType() == TRIANGLE_ELEMENT)
					{
						nids.push_back(nids.back());
					}

					Vector4d N;
					N(0) = (1 - s) * (1 - t) / 4.0;
					N(1) = (1 + s) * (1 - t) / 4.0;
					N(2) = (1 + s) * (1 + t) / 4.0;
					N(3) = (1 - s) * (1 + t) / 4.0;

					Vector4d X;
					Vector4d Y;
					Vector4d Z;
					for (int i = 0; i < 4; ++i)
					{
						const Vector3d& dis_local = T * ref_nodes[nids[i]].Displacement().block(0,0,3,1);
						X(i) = dis_local.x();
						Y(i) = dis_local.y();
						Z(i) = dis_local.z();
					}

					TDisplacement res_local;
					res_local.x() = N.transpose() * X;
					res_local.y() = N.transpose() * Y;
					res_local.z() = N.transpose() * Z;

					TDisplacement res_global;
					res_global.setZero();
					res_global.block(0,0,3,1) = T.transpose() * res_local.block(0, 0, 3, 1);

					return res_global;
				}
				void				UpdateLocalCoordinateSystem()
				{
					vector<int> nids = NodeIds();
					const TCoordinate& n1 = ref_nodes[nids[0]].Coordinate();
					const TCoordinate& n2 = ref_nodes[nids[1]].Coordinate();
					const TCoordinate& n3 = ref_nodes[nids[2]].Coordinate();

					Vector3d lx = n2 - n1;
					Vector3d ly = n3 - n2;
					Vector3d lz = Fork_Multi<Vector3d, Vector3d>(lx, ly);
					ly = Fork_Multi<Vector3d, Vector3d>(lz, lx);
					
					Normalizer<Vector3d>(lx);
					Normalizer<Vector3d>(ly);
					Normalizer<Vector3d>(lz);

					m_local_coord_system.block(0, 0, 1, 3) = lx.transpose();
					m_local_coord_system.block(1, 0, 1, 3) = ly.transpose();
					m_local_coord_system.block(2, 0, 1, 3) = lz.transpose();
				}
				Matrix3d&			LocalCoorSystem() { return m_local_coord_system; }
				const Matrix3d&		LocalCoorSystem() const { return m_local_coord_system; }

				MATRTIX_SINGLE_STIFFNESS&			SK() { return m_single_stiffness; }
				const MATRTIX_SINGLE_STIFFNESS&		SK() const { return m_single_stiffness; }

				TIntegrationPoint&			IP(int index) { return m_IP[index]; }
				const TIntegrationPoint&	IP(int index) const { return m_IP[index]; }
				void						UpdateIPCoordinate()
				{
					int nid1 = NodeId(0);
					int nid2 = NodeId(1);
					int nid3 = NodeId(2);
					int nid4 = NodeId(3);
					TNodeBase& node1 = ref_nodes[nid1];
					TNodeBase& node2 = ref_nodes[nid2];
					TNodeBase& node3 = ref_nodes[nid3];
					TNodeBase& node4 = ref_nodes[nid4];

					Eigen::VectorXd coord;
					coord.resize(12);
					coord.block(0, 0, 3, 1) = node1.Coordinate();
					coord.block(3, 0, 3, 1) = node2.Coordinate();
					coord.block(6, 0, 3, 1) = node3.Coordinate();
					coord.block(9, 0, 3, 1) = node4.Coordinate();

					//	更新单元积分点处的坐标信息
					for (int is = 0; is < IP_COUNT_2D; ++is)
					{
						//	计算积分点处的坐标时假定变形后的平面单元仍然保持平面形式
						Eigen::MatrixXd N = M_SF_RECTANGLE_PLANE(S_IP_2D[is],T_IP_2D[is]);
						IP(is).Coordinate() = N * coord;
					}
					//	最后一个积分点是单元的中心点
					IP(IP_COUNT_2D).Coordinate() = M_SF_RECTANGLE_PLANE(0, 0) * coord;
				}
			public:
				void				RefreshDatas()
				{
					UpdateLocalCoordinateSystem();
					vector<int> nids = NodeIds();
					if (MeshElementType() == TRIANGLE_ELEMENT)
					{
						nids.push_back(nids.back());
					}
					const Matrix3d& T = m_local_coord_system;

					const TCoordinate& n0 = ref_nodes[nids[0]].Coordinate();
					const TCoordinate& n1 = ref_nodes[nids[1]].Coordinate();
					const TCoordinate& n2 = ref_nodes[nids[2]].Coordinate();
					const TCoordinate& n3 = ref_nodes[nids[3]].Coordinate();

					m_side_a = Distance_2pt(n0, n1);
					m_side_b = Distance_2pt(n1, n2);
					m_side_length = (m_side_a + m_side_b) / 2;

					double a = Calculate_Area_Triangle(n0, n1, n2);
					double b = Calculate_Area_Triangle(n2, n3, n0);
					m_area = a + b;

					UpdateIPCoordinate();
				}
			private:
				vector<int*>		m_vec_nids;						//	Node Ids of this element
				int					m_part_id;						//	Part Id of this element
				int*				m_p_local_id;					//	Id of this element
				int					m_eid_global;					//	Global Id of this element
				ANALYSIS_ELEMENT_TYPE		m_elem_type;			//	Element type
			private:
				double				m_side_length;					//	Side Length of this element
				double				m_side_a;						//	Length of side A
				double				m_side_b;						//	Length of side B
				double				m_area;							//	Area of this element
				double				m_thickness;					//	Thickness of this element
			private:
				Matrix3d			m_local_coord_system;			//	Local coordinate system of the element
			private:
				MATRTIX_SINGLE_STIFFNESS		m_single_stiffness;				//	Single Stiffness Matrix of this element
				vector<TIntegrationPoint>	m_IP;					//	Integration Points of this element
			private:
				vector<TNodeBase>&	ref_nodes;						//	Node Data
			};

			//	FemMesh
			template<typename TNode, typename TElement>
			class TMeshCoreTemplate
			{
			public:
				TMeshCoreTemplate() {}
				~TMeshCoreTemplate()
				{
					for (vector<TNode>::iterator iterNode = m_nodes.begin();
						iterNode != m_nodes.end(); ++iterNode)
					{
						iterNode->Dispose();
					}
					m_nodes.clear();

					for (vector<TElement>::iterator iterElem = m_elements.begin();
						iterElem != m_elements.end(); ++iterElem)
					{
						iterElem->Dispose();
					}
					m_elements.clear();
				}
			public:
				void				Initialize()
				{
					m_nodes.clear();
					m_elements.clear();
				}
			public:
				void				InsertNode(const TCoordinate& pt, int global_nid = -1)
				{
					TNode node(pt, global_nid);
					m_nodes.push_back(node);
					m_nodes.back().Id() = (int)(m_nodes.size() - 1);
					m_set_nids.insert(m_nodes.back().Id());
				}
				void				DeleteNode(int nid)
				{
					//	节点是悬空点才可以删掉
					if (m_nodes[nid].AdjElementCount() == 0)
					{
						m_nodes.back().Id() = m_nodes[nid].Id();
						swap(m_nodes[nid], m_nodes.back());
						m_nodes.back().Dispose();
						m_nodes.pop_back();
					}
				}
				int					NodeCount() const {	return (int)(m_nodes.size());}
				TNode&				Node(int nid_local)
				{
					return m_nodes[nid_local];
				}
				const TNode&		Node(int nid_local) const
				{
					return m_nodes[nid_local];
				}
				const set<int>&		GetNodeIdsByAll() const
				{
					return m_set_nids;
				}
				set<int>			GetNodeIdsByPart(int pid)
				{
					set<int> res;
					res.clear();
					for (const TElement& element : m_elements)
					{
						if (element.PartId() == pid)
						{
							for (int nid : element.NodeIds())
							{
								res.insert(nid);
							}
						}
					}

					return res;
				}
		
				void				InsertElement(const vector<int>& nids, int id_global, int part_id)
				{
					m_elements.push_back(TElement(m_nodes));
					TElement& element = m_elements.back();
					element.Id() = (int)(m_elements.size() - 1);
					element.IdGlobal() = id_global;
					element.PartId() = part_id;
					for (int nid : nids)
					{						
						element.InsertNode(&(m_nodes[nid].Id()));

						TNode& node = Node(nid);
						node.InsertAdjElement(&(element.Id()));
					}
					element.RefreshDatas();
					m_set_eids.insert(element.Id());
				}
				void				DeleteElement(int eid)
				{
					TElement& elem = m_elements[eid];
					for (int i = 0; i < elem.NodeCount(); ++i)
					{
						TNode& node = m_nodes[elem.NodeId(i)];
						node.DeleteAdjElement(&(elem.Id()));
					}
					m_elements.back().Id() = m_elements[eid].Id();

					swap(m_elements[eid], m_elements.back());
					m_elements.back().Dispose();
					m_elements.pop_back();
				}
				int					ElementCount() const { return (int)(m_elements.size()); }
				TElement&			Element(int eid_local)
				{
					return m_elements[eid_local];
				}
				const TElement&		Element(int eid_local) const
				{
					return m_elements[eid_local];
				}
				void				AddSeparateElement(int eid)
				{
					TElement& element = Element(eid);
					//	只有是可断裂的PD区域，才能对网格进行PD网格的离散
					if (element.AnalysisElementType() == PD_ELEMENT)
					{					
						const vector<int>& nids_old = element.NodeIds();
						vector<int*> nids_new;
						nids_new.clear();
						for (int nid : nids_old)
						{
							int adjEleCount = m_nodes[nid].AdjElementCount();
							if (adjEleCount > 1)
							{							
								InsertNode(m_nodes[nid].Coordinate(), element.IdGlobal() * 100 + m_nodes[nid].IdGlobal());
							
								//	新节点添加单元拓扑信息
								m_nodes.back().InsertAdjElement(&element.Id());
								//	老节点与单元解除拓扑关系
								m_nodes[nid].DeleteAdjElement(&element.Id());

								//	把老节点的边界条件复制到新节点上去
								TNodeBase& node_new = m_nodes.back();
								node_new.BoundarySpcNode() = m_nodes[nid].BoundarySpcNode();
								node_new.InitVelocity() = m_nodes[nid].InitVelocity();
								node_new.BoundaryPreMotion() = m_nodes[nid].BoundaryPreMotion();
			
								node_new.LoadNodePoint() = m_nodes[nid].LoadNodePoint();
								node_new.ScaleLoadNodePoints(1.0 / adjEleCount);
								m_nodes[nid].ScaleLoadNodePoints(1.0 - 1.0 / adjEleCount);
								//	新节点编号放置于集合中用于单元信息更新
								nids_new.push_back(&(m_nodes.back().Id()));
								//	新节点编号放置到节点编号集合中
								m_set_nids.insert(node_new.Id());
							}
							else
							{
								nids_new.push_back(&(m_nodes[nid].Id()));
								//	新节点编号放置到节点编号集合中
								m_set_nids.insert(m_nodes[nid].Id());
							}							
						}
						element.UpdateNodes(nids_new);
					}
				}

				const set<int>&		GetElementIdsByAll() const
				{
					return m_set_eids;
				}
				set<int>			GetElementIdsByPart(int pid)
				{
					set<int> res;
					res.clear();
					for (const TElement& element : m_elements)
					{
						if (element.PartId() == pid)
						{
							res.insert(element.Id());
						}
					}

					return res;
				}
				set<int>			GetFemElementIdsByPart(int pid)
				{
					set<int> res;
					res.clear();
					for (const TElement& element : m_elements)
					{
						if ((element.PartId() == pid) || (element.AnalysisElementType() == FEM_ELEMENT))
						{
							res.insert(element.Id());
						}
					}

					return res;
				}

				set<int>			GetAdjElements(int eid, int iter_num = 1)
				{
					set<int> res;
					res.clear();
					res.insert(eid);
					const int pid = Element(eid).PartId();

					set<int> adj_node_ids;
					adj_node_ids.clear();
					while (iter_num >= 1)
					{
						for (int eid_find : res)
						{
							const TElement& element = Element(eid_find);
							vector<int> nids = element.NodeIds();
							for (int nid : nids)
							{
								adj_node_ids.insert(nid);
							}
						}

						for (int adj_nid : adj_node_ids)
						{
							const TNode& node = Node(adj_nid);
							set<int> adj_eids = node.AdjElementIds();

							for (int adj_eid : adj_eids)
							{
								//	如果找到的单元跟种子单元在同一个PART下，则视为相邻单元
								if (Element(adj_eid).PartId() == pid)
								{
									res.insert(adj_eid);
								}
							}
						}

						iter_num--;
					}

					return res;
				}
			private:
				vector<TNode>		m_nodes;
				vector<TElement>	m_elements;
			private:
				set<int>			m_set_nids;
				set<int>			m_set_eids;
			};

			typedef TMeshCoreTemplate<TNodeBase, TElementBase> TFemMeshCore;

			//	*MAT_ELASTIC/PLASTIC/
			class TMaterial
			{
			public:
				TMaterial(int mat_id = -1) : m_id(mat_id) {/*Do nothing*/ }
				~TMaterial() {/*Do nothing*/ }
			public:
				int& Id() { return m_id; }
				int						Id() const { return m_id; }
				string& Name() { return m_name; }
				string					Name() const { return m_name; }

				void					InsertMatData(string mat_name, double mat_value)
				{
					m_mat_datas.insert(pair<string, double>(mat_name, mat_value));
				}
				double					GetMatValue(string mat_name) const
				{
					map<string, double>::const_iterator iter = m_mat_datas.find(mat_name);
					if (iter != m_mat_datas.end())
					{
						return iter->second;
					}
					else
					{
						return 1E20;
					}
				}
			public:
				void					Initialize()
				{
					m_id = -1;
					m_mat_datas.clear();
				}
			private:
				map<string, double>		m_mat_datas;		//	<"ELASTIC", 200000>...
				int						m_id;				//	MAT ID
				string					m_name;				//	MAT NAME
			};

			//	*SECTION_SHELL/SOLID
			class TSection
			{
			public:
				TSection(int prop_id = -1) : m_id(prop_id) {/*Do nothing*/ }
				~TSection() {/*Do nothing*/ }
			public:
				int& Id() { return m_id; }
				int						Id() const { return m_id; }
				string& Type() { return m_type; }
				string					Type() const { return m_type; }

				void					InsertSectionData(string prop_name, double prop_value)
				{
					m_prop_datas.insert(pair<string, double>(prop_name, prop_value));
				}
				double					GetSectionValue(string prop_name) const
				{
					map<string, double>::const_iterator iter = m_prop_datas.find(prop_name);
					if (iter != m_prop_datas.end())
					{
						return iter->second;
					}
					else
					{
						return -10E10;
					}
				}
			public:
				void					Initialize()
				{
					m_id = -1;
					m_prop_datas.clear();
				}
			private:
				map<string, double>		m_prop_datas;	//	<"THICKNESS", 1.0>...
				int						m_id;			//	Section ID
				string					m_type;			//	SHELL or SOLID
			};  

			//	*PART
			class TPart
			{
			public:
				TPart(int part_id = -1, int mat_id = -1, int sec_id = -1)
				{
					m_pid = part_id;
					m_mid = mat_id;
					m_secid = sec_id;
				}
				~TPart()
				{
				}
			public:
				void				Initialize()
				{
					m_pid = -1;
					m_mid = -1;
					m_secid = -1;
					m_set_ids.clear();
				}
			public:
				int&				Id() { return m_pid; }
				int					Id() const { return m_pid; }
				int&				MaterialId() { return m_mid; }
				int					MaterialId() const { return m_mid; }
				int&				SectionId() { return m_secid; }
				int					SectionId() const { return m_secid; }
			public:
				const set<int>&		GetElementIds() const { return m_set_ids; }
				void				AddElementId(int nid) { m_set_ids.insert(nid); }
				void				AddElementId(const set<int>& eids)
				{
					for (const int& nid : eids)
					{
						m_set_ids.insert(nid);
					}
				}
				void				AddElementId(const vector<int>& eids)
				{
					for (const int& nid : eids)
					{
						m_set_ids.insert(nid);
					}
				}
			private:
				int					m_pid;			//	ID of part
				int					m_mid;			//	ID of material
				int					m_secid;		//	ID of section
				set<int>			m_set_ids;		//	Element of the part
			};

			//	*DEFINE_CURVE
			class TCurve
			{
			public:
				TCurve(int cur_id = -1) : m_id(cur_id) {/*Do nothing*/ }
				~TCurve() {/*Do nothing*/ }
			public:
				int& Id() { return m_id; }
				int			Id() const { return m_id; }

				double		GetValueByX(double x) const
				{
					double resY = 10E-10;
					if (!m_pts.empty())
					{
						resY = m_pts[0].y();
					}
					if (x < m_pts.front().x() + ERR_VALUE)
					{
						resY = m_pts.front().y();
					}
					else if (x > m_pts.back().x() - ERR_VALUE)
					{
						resY = m_pts.back().y();
					}
					else
					{
						for (int loop = 0; loop < (int)(m_pts.size()) - 1; ++loop)
						{
							double x_front = m_pts[loop].x() - ERR_VALUE;
							double x_back = m_pts[loop + 1].x() + ERR_VALUE;
							if ((x > x_front) && (x < x_back))
							{
								double x1 = m_pts[loop].x();
								double x2 = m_pts[loop + 1].x();
								double y1 = m_pts[loop].y();
								double y2 = m_pts[loop + 1].y();

								resY = y1 + (y2 - y1) / (x2 - x1) * (x - x1);
								break;
							}
						}
					}
					return resY;
				}
				void		AddPoint(double x, double y)
				{
					m_pts.push_back(Vector3d(x, y, 0));
				}
				double		GetXMax()
				{
					return m_pts.back().x();
				}
			public:
				void		Initialize()
				{
					m_id = -1;
					m_pts.clear();
				}
			private:
				vector<Vector3d>	m_pts;
				int					m_id;
			};
									
			//	*ELEMENT_SEATBELT
			class TCrevice
			{
			public:
				TCrevice(Vector3d p1 = Vector3d(), Vector3d p2 = Vector3d())
					: m_start(p1), m_end(p2)
				{ /*Do nothing*/
				}
				~TCrevice() { /*Do nothing*/ }
			public:
				Vector3d&		Start() { return m_start; }
				const Vector3d& Start() const { return m_start; }
				Vector3d&		End() { return m_end; }
				const Vector3d	End() const { return m_end; }
			private:
				Vector3d		m_start;
				Vector3d		m_end;
			};
		}
	}
}

#endif