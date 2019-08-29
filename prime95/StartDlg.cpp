// StartDlg.cpp : implementation file
//

#include "stdafx.h"
#include "Prime95.h"
#include "StartDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CStartDlg dialog


CStartDlg::CStartDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CStartDlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(CStartDlg)
	m_thread = 1;
	m_all_threads = 1;
	//}}AFX_DATA_INIT
}

void CStartDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CStartDlg)
	DDX_Check(pDX, IDC_RADIO1, m_all_threads);
	DDX_Control(pDX, IDC_THREAD_TEXT, c_thread_text);
	DDX_Control(pDX, IDC_THREAD, c_thread);
	DDX_Text(pDX, IDC_THREAD, m_thread);
	DDV_MinMaxUInt(pDX, m_thread, 1, WORKER_THREADS_ACTIVE > NUM_WORKER_THREADS ? WORKER_THREADS_ACTIVE : NUM_WORKER_THREADS);
	//}}AFX_DATA_MAP
	c_thread_text.EnableWindow (!m_all_threads);
	c_thread.EnableWindow (!m_all_threads);
}


BEGIN_MESSAGE_MAP(CStartDlg, CDialog)
	//{{AFX_MSG_MAP(CStartDlg)
	ON_BN_CLICKED(IDC_RADIO1, OnBnClickedRadio1)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()


/////////////////////////////////////////////////////////////////////////////
// CStartDlg message handlers

void CStartDlg::OnBnClickedRadio1()
{
	UpdateData ();
	c_thread_text.EnableWindow (!m_all_threads);
	c_thread.EnableWindow (!m_all_threads);
}

