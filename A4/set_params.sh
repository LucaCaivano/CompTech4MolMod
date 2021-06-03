#!/bin/sh

sed -e "s/%%DESTDIR%%/c2h6_${19}.chk/g" c2h6.com.template | \
sed -e "s/%%D1%%/$1/g" | \
sed -e "s/%%A1%%/$2/g" | \
sed -e "s/%%T1%%/$3/g" | \
sed -e "s/%%D2%%/$4/g" | \
sed -e "s/%%A2%%/$5/g" | \
sed -e "s/%%T2%%/$6/g" | \
sed -e "s/%%D3%%/$7/g" | \
sed -e "s/%%D4%%/$8/g" | \
sed -e "s/%%A3%%/$9/g" | \
sed -e "s/%%D5%%/${10}/g" | \
sed -e "s/%%A4%%/${11}/g" | \
sed -e "s/%%T3%%/${12}/g" | \
sed -e "s/%%D6%%/${13}/g" | \
sed -e "s/%%A5%%/${14}/g" | \
sed -e "s/%%T4%%/${15}/g" | \
sed -e "s/%%D7%%/${16}/g" | \
sed -e "s/%%A6%%/${17}/g" | \
sed -e "s/%%T5%%/${18}/g" > c2h6_${19}.com
