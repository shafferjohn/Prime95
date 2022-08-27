package main

import (
	"bytes"
	"fmt"
	"io/ioutil"
	"log"
	"net/http"
	"os"
	"strings"

	"github.com/PuerkitoBio/goquery"
)

type Producer struct {
	Rank         string
	TeamName     string
	TotalGHzDays string
	Day90        string
	Day30        string
	Day7         string
	Day1         string
	TF           string
	P_1          string
	LL_PRP       string
	DC           string
	ECM          string
	ECM_F        string
}

func (p Producer) String() string {
	return fmt.Sprintf("|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|", p.Rank, p.TeamName, p.TotalGHzDays, p.Day90, p.Day30, p.Day7, p.Day1, p.TF, p.P_1, p.LL_PRP, p.DC, p.ECM, p.ECM_F)
}

func main() {
	res, err := http.Get("https://www.mersenne.org/report_top_teams/")
	if err != nil {
		log.Fatal(err)
	}
	defer res.Body.Close()
	if res.StatusCode != 200 {
		log.Fatal("status code error: %d %s", res.StatusCode, res.Status)
	}

	doc, err := goquery.NewDocumentFromReader(res.Body)
	if err != nil {
		log.Fatal(err)
	}

	ths := []string{}
	data := make([]Producer, 0)
	doc.Find("thead tr:last-child th").Each(func(i int, s *goquery.Selection) {
		ths = append(ths, s.Text())
	})
	caption, _ := doc.Find("caption").Html()
	caption = strings.ReplaceAll(caption, "<br>", "\n\n")
	doc.Find("tbody").Last().Find("tr").Each(func(i int, s *goquery.Selection) {
		p := Producer{}
		s.Find("td").Each(func(j int, ss *goquery.Selection) {
			arrow := ss.Find("img").AttrOr("src", "")
			if strings.Contains(arrow, "up") {
				arrow = "⬆ "
			} else if strings.Contains(arrow, "dn") {
				arrow = "⬇ "
			}

			switch j {
			case 0:
				p.Rank = ss.Text()
			case 1:
				p.TeamName = ss.Text()
				link := ss.Find("a:first-child").AttrOr("href", "")
				if link != "" {
					p.TeamName = fmt.Sprintf("[%s](%s)", p.TeamName, link)
				}
			case 2:
				p.TotalGHzDays = ss.Text()
			case 3:
				p.Day90 = fmt.Sprintf("%v%v", arrow, ss.Text())
			case 4:
				p.Day30 = fmt.Sprintf("%v%v", arrow, ss.Text())
			case 5:
				p.Day7 = fmt.Sprintf("%v%v", arrow, ss.Text())
			case 6:
				p.Day1 = fmt.Sprintf("%v%v", arrow, ss.Text())
			case 7:
				p.TF = ss.Text()
			case 8:
				p.P_1 = ss.Text()
			case 9:
				p.LL_PRP = ss.Text()
			case 10:
				p.DC = ss.Text()
			case 11:
				p.ECM = ss.Text()
			case 12:
				p.ECM_F = ss.Text()
			}
		})
		data = append(data, p)
	})

	content := "\n" + caption + "\n\n"
	content += "|" + strings.Join(ths, "|") + "|\n"
	for i := 0; i < len(ths); i++ {
		content += "|----"
	}
	content += "|\n"
	for _, p := range data {
		content += p.String() + "\n"
	}

	filename := "README.md"
	md, err := ioutil.ReadFile(filename)
	if err != nil {
		log.Fatal(err)
	}
	start := []byte("<!-- PrimeNet Top Rank start -->")
	before := md[:bytes.Index(md, start)+len(start)]
	end := []byte("<!-- PrimeNet Top Rank end -->")
	after := md[bytes.Index(md, end):]

	newMD := bytes.NewBuffer(nil)
	newMD.Write(before)
	newMD.WriteString("\n" + content + "\n")
	newMD.Write(after)
	err = ioutil.WriteFile(filename, newMD.Bytes(), os.ModeAppend)
	if err != nil {
		log.Fatal(err)
	}
}
