# SIPRI raw data


How to download the SIPRI data:

```text
curl http://armstrade.sipri.org/armstrade/html/export_trade_register.php --compressed \
    --data 'low_year=2014' \
    --data 'high_year=2014' \
    --data 'seller_country_code=' \
    --data 'buyer_country_code=' \
    --data 'armament_category_id=any' \
    --data 'buyers_or_sellers=sellers' \
    --data 'filetype=csv' \
    --data 'include_open_deals=on' \
    --data 'sum_deliveries=on' \
    --data 'Submit4=Download' \
> sipri-arms-by-seller-2014.csv
```
