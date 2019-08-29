// EcmDlg.cpp : implementation file
//
//  Copyright 2000-2019 Mersenne Research, Inc. All rights reserved.
//

#include "stdafx.h"
#include "Prime95.h"
#include "EcmDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CEcmDlg dialog


CEcmDlg::CEcmDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CEcmDlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(CEcmDlg)
	m_thread = 1;
	m_k = 1.0;
	m_b = 2;
	m_n = 1277;
	m_c = -1;
	m_bound1 = 850.0e6;
	m_bound2 = 0.0;
	m_num_curves = 10;
	//}}AFX_DATA_INIT
}


void CEcmDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CEcmDlg)
	DDX_Control(pDX, IDC_THREAD_TEXT, c_thread_text);
	DDX_Control(pDX, IDC_THREAD, c_thread);
	DDX_Text(pDX, IDC_THREAD, m_thread);
	DDV_MinMaxUInt(pDX, m_thread, 1, NUM_WORKER_THREADS);
	DDX_Text(pDX, IDC_P4, m_k);
	DDX_Text(pDX, IDC_P1, m_b);
	DDV_MinMaxUInt (pDX, m_b, 2, 1000000000);
	DDX_Text(pDX, IDC_P5, m_n);
	DDV_MinMaxUInt (pDX, m_n, 1,
			CPU_FLAGS & CPU_SSE2 ? MAX_PRIME_SSE2 : MAX_PRIME);
	DDX_Text(pDX, IDC_P6, m_c);
	DDX_Text(pDX, IDC_P2, m_bound1);
	DDX_Text(pDX, IDC_P3, m_bound2);
	DDX_Text(pDX, IDC_NUM_CURVES, m_num_curves);
	DDV_MinMaxUInt(pDX, m_num_curves, 1, 1000000);
	//}}AFX_DATA_MAP
	c_thread_text.EnableWindow (NUM_WORKER_THREADS > 1);
	c_thread.EnableWindow (NUM_WORKER_THREADS > 1);
}


BEGIN_MESSAGE_MAP(CEcmDlg, CDialog)
	//{{AFX_MSG_MAP(CEcmDlg)
		// NOTE: the ClassWizard will add message map macros here
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CEcmDlg message handlers
