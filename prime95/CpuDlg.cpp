// CpuDlg.cpp : implementation file
//

#include "stdafx.h"
#include "Prime95.h"
#include "CpuDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CCpuDlg dialog


CCpuDlg::CCpuDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CCpuDlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(CCpuDlg)
	m_hours = 0;
	m_start_time = _T("");
	m_end_time = _T("");
	m_day_memory = 0;
	m_night_memory = 0;
	m_cpu_info = _T("");
	//}}AFX_DATA_INIT
}


void CCpuDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CCpuDlg)
	DDX_Text(pDX, IDC_HOURS, m_hours);
	DDV_MinMaxUInt(pDX, m_hours, 1, 24);
	DDX_Text(pDX, IDC_START_TIME, m_start_time);
	DDX_Text(pDX, IDC_END_TIME, m_end_time);
	DDX_Text(pDX, IDC_DAY_MEMORY, m_day_memory);
	DDV_MinMaxUInt(pDX, m_day_memory, 8, (UINT) (0.9 * physical_memory ()));
	DDX_Text(pDX, IDC_NIGHT_MEMORY, m_night_memory);
	DDV_MinMaxUInt(pDX, m_night_memory, 8, (UINT) (0.9 * physical_memory ()));
	DDX_Text(pDX, IDC_CPU_INFO, m_cpu_info);

	DDX_Control(pDX, IDC_DAY_MEMORY_TEXT, c_day_memory_text);
	DDX_Control(pDX, IDC_DAY_MEMORY, c_day_memory);
	DDX_Control(pDX, IDC_NIGHT_MEMORY_TEXT, c_night_memory_text);
	DDX_Control(pDX, IDC_NIGHT_MEMORY, c_night_memory);
	DDX_Control(pDX, IDC_START_TIME_TEXT, c_start_time_text);
	DDX_Control(pDX, IDC_START_TIME, c_start_time);
	DDX_Control(pDX, IDC_END_TIME_TEXT, c_end_time_text);
	DDX_Control(pDX, IDC_END_TIME, c_end_time);
	//}}AFX_DATA_MAP
	c_day_memory_text.EnableWindow (m_memory_editable);
	c_day_memory.EnableWindow (m_memory_editable);
	c_night_memory_text.EnableWindow (m_memory_editable);
	c_night_memory.EnableWindow (m_memory_editable);
	c_start_time_text.EnableWindow (m_memory_editable);
	c_start_time.EnableWindow (m_memory_editable);
	c_end_time_text.EnableWindow (m_memory_editable);
	c_end_time.EnableWindow (m_memory_editable);
}


BEGIN_MESSAGE_MAP(CCpuDlg, CDialog)
	//{{AFX_MSG_MAP(CCpuDlg)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CCpuDlg message handlers

