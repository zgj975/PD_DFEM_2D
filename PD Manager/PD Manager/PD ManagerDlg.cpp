
// PD ManagerDlg.cpp : 实现文件
//

#include "stdafx.h"
#include "PD Manager.h"
#include "PD ManagerDlg.h"
#include "afxdialogex.h"
#include <string>
#include "peridynamic/pd_product.hxx"
using namespace DLUT::SAE::PERIDYNAMIC;
using namespace std;

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// 用于应用程序“关于”菜单项的 CAboutDlg 对话框

class CAboutDlg : public CDialogEx
{
public:
	CAboutDlg();

// 对话框数据
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_ABOUTBOX };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

// 实现
protected:
	DECLARE_MESSAGE_MAP()	
};

CAboutDlg::CAboutDlg() : CDialogEx(IDD_ABOUTBOX)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialogEx)
END_MESSAGE_MAP()

// CPDManagerDlg 对话框
CPDManagerDlg::CPDManagerDlg(CWnd* pParent /*=NULL*/)
	: CDialogEx(IDD_PDMANAGER_DIALOG, pParent)
	, m_fp_information(_T(""))
	, m_ratio_of_delta_meshsize(3.0)
	, m_time_step(1E-7)
	, m_plot_frames(100)
	, m_radio_analysis_type(0)
	, m_iterator_nums(1000)
	, m_load_step(1)
	, m_radio_horizon_type(0)
	, m_d_const_horizon(3.0)
	, m_max_iter_nums(100)
	, m_convergence_factor(0.001)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CPDManagerDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_EDIT_FILE_PD_INFORMATION, m_fp_information);
	DDX_Text(pDX, IDC_EDIT_RATIO_D_M, m_ratio_of_delta_meshsize);
	DDX_Text(pDX, IDC_EDIT_TIME_STEP, m_time_step);
	DDX_Text(pDX, IDC_EDIT_PLOT_FRAME, m_plot_frames);
	DDX_Radio(pDX, IDC_RADIO_IMPLICIT, m_radio_analysis_type);
	DDX_Text(pDX, IDC_EDIT_ITERATOR_NUMS, m_iterator_nums);
	DDX_Text(pDX, IDC_EDIT_LOAD_STEP, m_load_step);
	DDX_Radio(pDX, IDC_RADIO_CONSTANT_HORIZON, m_radio_horizon_type);
	DDX_Text(pDX, IDC_EDIT_CONSTANT_HORIZON, m_d_const_horizon);
	DDX_Text(pDX, IDC_EDIT_MAX_ITERATOR_NUMS, m_max_iter_nums);
	DDX_Text(pDX, IDC_EDIT_CONVERGENCE_FACTOR, m_convergence_factor);
}

BEGIN_MESSAGE_MAP(CPDManagerDlg, CDialogEx)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_BUTTON_PD_INFORMATION, &CPDManagerDlg::OnBnClickedButtonPdInformation)
	ON_BN_CLICKED(IDC_RADIO_IMPLICIT, &CPDManagerDlg::OnBnClickedRadioImplicit)
	ON_BN_CLICKED(IDOK, &CPDManagerDlg::OnBnClickedOk)
	ON_BN_CLICKED(IDC_RADIO_EXPLICIT, &CPDManagerDlg::OnBnClickedRadioExplicit)
	ON_BN_CLICKED(IDC_RADIO_CONSTANT_HORIZON, &CPDManagerDlg::OnBnClickedRadioConstantHorizon)
	ON_BN_CLICKED(IDC_RADIO_CONSTANT_M, &CPDManagerDlg::OnBnClickedRadioConstantM)
END_MESSAGE_MAP()


// CPDManagerDlg 消息处理程序

BOOL CPDManagerDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// 将“关于...”菜单项添加到系统菜单中。

	// IDM_ABOUTBOX 必须在系统命令范围内。
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		BOOL bNameValid;
		CString strAboutMenu;
		bNameValid = strAboutMenu.LoadString(IDS_ABOUTBOX);
		ASSERT(bNameValid);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// 设置此对话框的图标。  当应用程序主窗口不是对话框时，框架将自动
	//  执行此操作
	SetIcon(m_hIcon, TRUE);			// 设置大图标
	SetIcon(m_hIcon, FALSE);		// 设置小图标

	// TODO: 在此添加额外的初始化代码

	return TRUE;  // 除非将焦点设置到控件，否则返回 TRUE
}

void CPDManagerDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialogEx::OnSysCommand(nID, lParam);
	}
}

// 如果向对话框添加最小化按钮，则需要下面的代码
//  来绘制该图标。  对于使用文档/视图模型的 MFC 应用程序，
//  这将由框架自动完成。

void CPDManagerDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // 用于绘制的设备上下文

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// 使图标在工作区矩形中居中
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// 绘制图标
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

//当用户拖动最小化窗口时系统调用此函数取得光标
//显示。
HCURSOR CPDManagerDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}

void CPDManagerDlg::OnBnClickedButtonPdInformation()
{
	// TODO: 在此添加控件通知处理程序代码
	BOOL isOpen = TRUE;
	CString filter = L"PD Information (*.key; *.k)|*.key;*.k|All Files (*.*)|*.*|";
	CFileDialog openFileDlg(isOpen, NULL, NULL, OFN_HIDEREADONLY | OFN_READONLY, filter, NULL);
	INT_PTR result = openFileDlg.DoModal();

	if (result == IDOK)
	{
		CString filePath = openFileDlg.GetPathName();
		CWnd::SetDlgItemTextW(IDC_EDIT_FILE_PD_INFORMATION, filePath);
	}
}

void CPDManagerDlg::OnBnClickedRadioImplicit()
{
	// TODO: 在此添加控件通知处理程序代码
	UpdateData(TRUE);
	
	GetDlgItem(IDC_EDIT_LOAD_STEP)->EnableWindow(TRUE);
	GetDlgItem(IDC_EDIT_MAX_ITERATOR_NUMS)->EnableWindow(TRUE);
	GetDlgItem(IDC_EDIT_CONVERGENCE_FACTOR)->EnableWindow(TRUE);

	GetDlgItem(IDC_EDIT_TIME_STEP)->EnableWindow(FALSE);
	GetDlgItem(IDC_EDIT_ITERATOR_NUMS)->EnableWindow(FALSE);
	GetDlgItem(IDC_EDIT_PLOT_FRAME)->EnableWindow(FALSE);
		
	UpdateData(FALSE);
}

void CPDManagerDlg::OnBnClickedRadioExplicit()
{
	// TODO: 在此添加控件通知处理程序代码
	UpdateData(TRUE);

	GetDlgItem(IDC_EDIT_LOAD_STEP)->EnableWindow(FALSE);
	GetDlgItem(IDC_EDIT_MAX_ITERATOR_NUMS)->EnableWindow(FALSE);
	GetDlgItem(IDC_EDIT_CONVERGENCE_FACTOR)->EnableWindow(FALSE);

	GetDlgItem(IDC_EDIT_TIME_STEP)->EnableWindow(TRUE);
	GetDlgItem(IDC_EDIT_ITERATOR_NUMS)->EnableWindow(TRUE);
	GetDlgItem(IDC_EDIT_PLOT_FRAME)->EnableWindow(TRUE);

	UpdateData(FALSE);
}

void CPDManagerDlg::OnBnClickedRadioConstantHorizon()
{
	// TODO: 在此添加控件通知处理程序代码
	UpdateData(TRUE);

	GetDlgItem(IDC_EDIT_CONSTANT_HORIZON)->EnableWindow(TRUE);
	GetDlgItem(IDC_EDIT_RATIO_D_M)->EnableWindow(FALSE);

	UpdateData(FALSE);
}

void CPDManagerDlg::OnBnClickedRadioConstantM()
{
	// TODO: 在此添加控件通知处理程序代码
	UpdateData(TRUE);

	GetDlgItem(IDC_EDIT_CONSTANT_HORIZON)->EnableWindow(FALSE);
	GetDlgItem(IDC_EDIT_RATIO_D_M)->EnableWindow(TRUE);

	UpdateData(FALSE);
}

void CPDManagerDlg::OnBnClickedOk()
{
	// TODO: 在此添加控件通知处理程序代码
	UpdateData(TRUE);
		
	DWORD dwAttrib = GetFileAttributes(m_fp_information);
	if (INVALID_FILE_ATTRIBUTES == dwAttrib || 0 != (dwAttrib & FILE_ATTRIBUTE_DIRECTORY))
	{
		MessageBox(L"You should select a valid LSDYNA file, first!", L"INFORMATION", MB_OK);
		return;
	}
	
	AllocConsole();
	FILE* fp_out;
	FILE* fp_in;
	freopen_s(&fp_out, "CONOUT$", "w+t", stdout);
	freopen_s(&fp_in, "CONIN$", "w+t", stdin);

	string fp_lsdyna = CW2A(m_fp_information.GetString());

	PdProduct<TSolverBeamPD> product;

	if (m_radio_horizon_type == 0)
	{
		product.SetUsedConstantHorizon(true);
		product.SetConstantHorizon(m_d_const_horizon);
	}
	else
	{
		product.SetUsedConstantHorizon(false);
		product.SetRatioOfHorizonMeshsize(m_ratio_of_delta_meshsize);
	}
	product.SetFilePath(fp_lsdyna);

	if (m_radio_analysis_type == 0)
	{
		product.SetLoadStep(m_load_step);
		product.SetConvergenceFactor(m_convergence_factor);
		product.SetMaxIteratorNums(m_max_iter_nums);

		product.ImplicitAnalysis();
	}
	else if (m_radio_analysis_type == 1)
	{
		product.SetTimeStep(m_time_step);
		product.SetIteratorNums(m_iterator_nums);
		product.SetPlotFrames(m_plot_frames);

		product.ExplicitAnalysis();
	}
	fclose(fp_in);
	fclose(fp_out);

	cout << endl << endl;
	system("pause");

	FreeConsole();	
}
