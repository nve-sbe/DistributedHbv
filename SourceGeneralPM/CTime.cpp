// $Id: ctime.C 399 2006-02-15 12:48:48Z bjg $
//      File: time.c

// Last edited by Trond Reitan, 20/4-2004
// Replaced form() with sprintf in operator<<

#include "CTime.h"
#include <string.h>
#include <fstream>
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>

int legalTime(const short h, const short m, const short s)
{
    return (h < 0 || h > 23 || m < 0 || m > 59 || s < 0 || s > 59) ?
           0 : 1;
}


extern void skipDelim(istream&);

totalSecondsType CTime::time2totSec(hourType h, minuteType m, secondType s) const
{
    /* Convert an hhmmss time (hour:minute:sec) to seconds.
    This function is not valid for time whose total seconds is greater
    than 2**sizeof(unsigned long). */
    if (!legalTime(h, m, s))
    {
        h = 24;
        m = 60;
        s = 60;
    };
    return (s + m * 60 + h * 60 * 60);

} /* time2totSec */


void CTime::totSec2time(hourType &h, minuteType &m, secondType &s) const
{
    /* Convert seconds to hhmmss time (hour:minute:sec). This function is
    not valid for time whose total seconds is greter
    than 2**sizeof(unsigned long). */

    h = tsec / 3600;
    m = (tsec % 3600) / 60;
    s = (tsec % 3600) % 60;

} /* totSec2time */

void CTime::now()
{
    time_t clk = time(0);
    tm* now = NULL;
    now = localtime(&clk);
    tsec = time2totSec(now->tm_hour, now->tm_min, now->tm_sec);
}

CTime::CTime(const CTime &t)
{
    tsec = t.tsec;
}


CTime::CTime(hourType h, minuteType m, secondType s)
{
    tsec = time2totSec(h, m, s);
}

CTime::CTime(const char *time)
{
    static char tem[3];
    static hourType hh = 0;
    static minuteType mm = 0;
    static secondType ss = 0;
    tem[2] = '\0';
    hh = (hourType)atoi(strncpy(tem, time, 2));
    mm = (minuteType)atoi(strncpy(tem, time + 2, 2));
    ss = 0;
    if (strlen(time) > 4)
    {
        ss = atoi(strncpy(tem, time + 4, 2));
    }
    if (legalTime(hh, mm, ss))
    {
        tsec = time2totSec(hh, mm, ss);
    }
    else
    {
        tsec = time2totSec(24, 60, 60);
    }
}

int CTime::legal() const
{
    static hourType hh = 0;
    static minuteType mm = 0;
    static secondType ss = 0;
    totSec2time(hh, mm, ss);
    return (::legalTime(hh, mm, ss));
}

totalSecondsType CTime::parseTime(istream& is) const
{
    /* parse time of the format hh:mm:ss
    Delimiter need not be ':',
    it can be any of the alphanumeric characters */

    unsigned h = 0, m = 0, s = 0;

    if (is.good())
    {
        skipDelim(is);
        is >> h;
        skipDelim(is);
        if (is.eof())
        {
            return 0;
        }

        is >> m;
        skipDelim(is);
        if (is.eof())
        {
            return 0;
        }
        is >> s;
    }
    if (!is.good())
    {
        return 0;
    }
    return CTime((hourType)h, m, s).tsec;
}

// end parseTime

CTime::CTime(istream &is)
{
    tsec = parseTime(is);
}


ostream& operator<<(ostream& os, const CTime &t)
{
    static hourType hh = 0;
    static minuteType mm = 0;
    static secondType ss = 0;
    static char buffer[10];

    t.totSec2time(hh, mm, ss);
    sprintf(buffer, "%02d:%02d", hh, mm);
    os << buffer;

    return os;
}

char *CTime::getChTime() const
{
    static char ctime[9];
    static hourType h = 0;
    static minuteType mn = 0;
    static secondType s = 0;

    totSec2time(h, mn, s);
    sprintf(ctime, "%02d:%02d", h, mn);
    return ctime;
}

hourType CTime::getHour() const
{
    static hourType hh = 0;
    static minuteType mm = 0;
    static secondType ss = 0;
    totSec2time(hh, mm, ss);
    return hh;
}

minuteType CTime::getMinute() const
{
    static hourType hh = 0;
    static minuteType mm = 0;
    static secondType ss = 0;
    totSec2time(hh, mm, ss);
    return mm;
}

secondType CTime::getSecond() const
{
    static hourType hh = 0;
    static minuteType mm = 0;
    static secondType ss = 0;
    totSec2time(hh, mm, ss);
    return ss;
}


int CTime::minofday() const
{
    return tsec / 60;
}
