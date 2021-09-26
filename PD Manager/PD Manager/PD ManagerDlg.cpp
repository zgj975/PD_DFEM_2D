
// PD ManagerDlg.cpp : ʵ���ļ�
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


// ����Ӧ�ó��򡰹��ڡ��˵���� CAboutDlg �Ի���

class CAboutDlg : public CDialogEx
{
public:
	CAboutDlg();

// �Ի�������
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_ABOUTBOX };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV ֧��

// ʵ��
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

// CPDManagerDlg �Ի���
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


// CPDManagerDlg ��Ϣ�������

BOOL CPDManagerDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// ��������...���˵�����ӵ�ϵͳ�˵��С�

	// IDM_ABOUTBOX ������ϵͳ���Χ�ڡ�
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

	// ���ô˶Ի����ͼ�ꡣ  ��Ӧ�ó��������ڲ��ǶԻ���ʱ����ܽ��Զ�
	//  ִ�д˲���
	SetIcon(m_hIcon, TRUE);			// ���ô�ͼ��
	SetIcon(m_hIcon, FALSE);		// ����Сͼ��

	// TODO: �ڴ���Ӷ���ĳ�ʼ������

	return TRUE;  // ���ǽ��������õ��ؼ������򷵻� TRUE
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

// �����Ի��������С����ť������Ҫ����Ĵ���
//  �����Ƹ�ͼ�ꡣ  ����ʹ���ĵ�/��ͼģ�͵� MFC Ӧ�ó���
//  �⽫�ɿ���Զ���ɡ�

void CPDManagerDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // ���ڻ��Ƶ��豸������

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// ʹͼ���ڹ����������о���
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// ����ͼ��
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

//���û��϶���С������ʱϵͳ���ô˺���ȡ�ù��
//��ʾ��
HCURSOR CPDManagerDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}

void CPDManagerDlg::OnBnClickedButtonPdInformation()
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
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
	// TODO: �ڴ���ӿؼ�֪ͨ����������
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
	// TODO: �ڴ���ӿؼ�֪ͨ����������
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
	// TODO: �ڴ���ӿؼ�֪ͨ����������
	UpdateData(TRUE);

	GetDlgItem(IDC_EDIT_CONSTANT_HORIZON)->EnableWindow(TRUE);
	GetDlgItem(IDC_EDIT_RATIO_D_M)->EnableWindow(FALSE);

	UpdateData(FALSE);
}

void CPDManagerDlg::OnBnClickedRadioConstantM()
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
	UpdateData(TRUE);

	GetDlgItem(IDC_EDIT_CONSTANT_HORIZON)->EnableWindow(FALSE);
	GetDlgItem(IDC_EDIT_RATIO_D_M)->EnableWindow(TRUE);

	UpdateData(FALSE);
}

void CPDManagerDlg::OnBnClickedOk()
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
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
